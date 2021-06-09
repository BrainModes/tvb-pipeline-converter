from Crypto.PublicKey import RSA
from Crypto.Cipher import PKCS1_OAEP
import pyAesCrypt
import sys

private_key_compute = sys.argv[1]
private_key_compute_pwd = sys.argv[2]
private_key_input_pwd = sys.argv[3]
private_key_input = sys.argv[4]
encrypted_password = sys.argv[5]
encrypted_input_data = sys.argv[6]

# Step 1: Decrypt private key private_key_compute.bin
encoded_key = open(private_key_compute, "rb").read()
private_key_compute = RSA.import_key(encoded_key, passphrase=private_key_compute_pwd)

# Step 2: Decrypt password private_key_input_pwd with decrypted private_key_compute
cipher_rsa = PKCS1_OAEP.new(private_key_compute)
private_key_input_pwd = cipher_rsa.decrypt(open(private_key_input_pwd, "rb").read())

# Step 3: Decrypt private key private_key_input.bin with private_key_input_pwd
encoded_key = open(private_key_input, "rb").read()
private_key_input = RSA.import_key(encoded_key, passphrase=private_key_input_pwd)

# Step 4: Decrypt password encrypted_password.bin for encrypted input data using private_key_input
cipher_rsa = PKCS1_OAEP.new(private_key_input)
session_key = cipher_rsa.decrypt(open(encrypted_password).read())

# Step 5: Decrypt the encrypted data with the AES session key
bufferSize = 64 * 1024
pyAesCrypt.decryptFile(encrypted_input_data, "/data/data.zip", session_key, bufferSize)

