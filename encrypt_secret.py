from Crypto.PublicKey import RSA
from Crypto.Cipher import PKCS1_OAEP
import sys

password = sys.argv[1]
output_dir = sys.argv[2]
public_key_file = sys.argv[3]


# Step 1: Load public RSA key generated on compute node
recipient_key = RSA.import_key(open(public_key_file).read())


# Step 2: Encrypt the password with the public RSA key
cipher_rsa = PKCS1_OAEP.new(recipient_key)
enc_password = cipher_rsa.encrypt(str.encode(password))
password = "" # overwrite password clear text in memory


# Step 3: Write out encrypted password
file_password = open(output_dir + "/" + "private_key_input_pwd.bin", "wb")
file_password.write(enc_password) 
file_password.close()
