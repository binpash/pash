#!/bin/bash

rsync -avz -e "ssh " /path/to/yourfile user@backupserver.com:/backup/
echo "backup for $(date) " | mail -s "backup complete" user@youremail.com
