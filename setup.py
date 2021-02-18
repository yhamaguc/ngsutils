import sys
import os
import re
import ftplib

from setuptools import setup


setup()

os.makedirs("assets", exist_ok=True)

ftp_server = "ftp.ebi.ac.uk"
ftp_path = "pub/databases/gencode/Gencode_human/latest_release"
ftp = ftplib.FTP(ftp_server)

ftp.login()
ftp.cwd(ftp_path)
ftp.nlst(".")

pattern = r"gencode.v([0-9])+.annotation.gtf.gz"
gtf = [f for f in ftp.nlst(".") if re.match(pattern, f)][0]

out_path = 'assets/gencode.annotation.gtf.gz'
sys.stderr.write(f"Now downloading... {ftp_server}/{ftp_path}/{gtf}")

with open(out_path, "wb") as f:
    ftp.retrbinary(f"RETR {gtf}", f.write)

os.system(f"gtf2sqlite {out_path}")
