Note: You need to have a recent copy of java to use the script.

Step 1: Edit config.cfg
-------

If you are running the script manually:
Add zID after
user=

leave password line blank:
password=
This will force the script to ask for your zPass when you run the script. It is not recommended to save your zPass here as it will be visible to anyone who has access to the file.

If you wish to automate the script please contact the servicedesk requesting a Token be created for your username and the RDMP number the token needs access to. See the config.cfg for where to put this token when you get it.

In the directory containing the script files run the following command to upload data to the Archive:

For Windows:
upload.bat "D:\path\to\your\local\directory" "/UNSW_RDS/D0000000/your/collection/name"
For Linux/Mac OS:
./upload.sh "/path/to/your/local/directory" "/UNSW_RDS/D0000000/your/collection/name"

To restore data to make it available for download:

For Windows:
restore.bat "/UNSW_RDS/D0000000/your/collection/name"
For Linux/Mac OS:
./restore.sh "/UNSW_RDS/D0000000/your/collection/name"

To download data from the archive:

For Windows:
download.bat "/UNSW_RDS/D0000000/your/collection/name" "D:\path\to\your\local\directory"
For Linux/Mac OS:
./download.sh "/UNSW_RDS/D0000000/your/collection/name" "/path/to/your/local/directory"

To download a single file:

For Windows:
download-file.bat "/UNSW_RDS/D0000000/your/collection/name/your_file.txt" "D:\path\to\your\local\directory\your_file.txt"
For Linux/Mac OS:
./download-file.sh "/UNSW_RDS/D0000000/your/collection/name/your_file.txt" "/path/to/your/local/directory/your_file.txt"

For Checksum Generation:
=========

Complete the above setup in config.cfg, then run the script as follows:

For Linux/Mac:
./getchecksums.sh "/UNSW_RDS/D000000X/path/to/files"

For Windows:
getchecksums.bat "/UNSW_RDS/D000000X/path/to/files"

This will create a file "checksums.csv" in the working directory containing CRC32 checksums for files in the specified folder/path in the Data Archive.
