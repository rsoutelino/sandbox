import imaplib
import email


imap_ssl_host = 'imap.gmail.com'
imap_ssl_port = 993
username = "someone@gmail.com"
password = "SomeOnePassword"
server = imaplib.IMAP4_SSL(imap_ssl_host, imap_ssl_port)

server.login(username, password)
server.select('INBOX')

status, data = server.search(None, 'ALL')

emails = []

for num in data[0].split():
    print(num)
    status, obj = server.fetch(num, '(RFC822)')
    emails.append(email.message_from_bytes(obj[0][1]))

for mail in emails:
    print(mail['From'])
    print(mail.as_string().split('\n\n')[-1])



