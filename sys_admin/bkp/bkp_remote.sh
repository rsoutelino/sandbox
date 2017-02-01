
rsync -Cravpe ssh /home/rsoutelino/ golfo:/bkp/rsoutelino/  --exclude-from '/home/rsoutelino/admin/bkp/bkp_exclude.txt'

# Every first business day of the month, please run the two 
#   scripts below. The first one will only list the files that
#   will be copied or deleted, and the second actually executes 
#   the task. Run the first, check things out and than run the 
#   second. 
 
# rsync -Cravpen ssh --delete /home/rsoutelino/ golfo:/bkp/rsoutelino/
# rsync -Cravpe  ssh --delete /home/rsoutelino/ golfo:/bkp/rsoutelino/
