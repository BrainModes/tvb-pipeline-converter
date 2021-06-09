from Crypto.PublicKey import RSA
from Crypto.Cipher import PKCS1_OAEP
import pyAesCrypt
import sys

public_key_results = sys.argv[1]
results_targz = sys.argv[2]
results_targz_aes = sys.argv[3]
encrypted_results_password = sys.argv[4]

# Step 1: create random password
pass_size = 40 # length of the password
password = ''.join(random.choice(string.ascii_letters + string.digits) for i in range(pass_size))


# Step 2: encrypt data using perviously generated password
bufferSize = 64 * 1024
pyAesCrypt.encryptFile(results_targz, results_targz_aes, password, bufferSize)


# Step 3: load public RSA key generated on data controller's computer
recipient_key = RSA.import_key(open(public_key_results).read())


# Step 4: encrypt the password with the public RSA key
cipher_rsa = PKCS1_OAEP.new(recipient_key)
enc_password = cipher_rsa.encrypt(str.encode(password))
password = "" # overwrite password clear text in memory


# Step 5: write out encrypted password
file_password = open(encrypted_results_password, "wb")
file_password.write(enc_password) 
file_password.close()

