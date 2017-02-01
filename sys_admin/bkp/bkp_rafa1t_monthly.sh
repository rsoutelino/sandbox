
# Every first business day of the month, please run the two 
#   scripts below. The first one will only list the files that
#   will be copied or deleted, and the second actually executes 
#   the task. Run the first, check things out and than run the 
#   second. 
 
# rsync -Cravpztun --delete /home/rsoutelino/ /media/RAFA500G/rsoutelino/
rsync -Cravpztu  --delete /home/rsoutelino/ /media/rsoutelino/RAFA1T/rsoutelino/
