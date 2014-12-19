#!/usr/bin/env python
import time
import sys
import smtplib
import socket



def set_message(info, imp):
    hostname = socket.gethostname()
    message = """Subject: %s: %s %s 
    """ %(hostname.upper(), imp.upper(), info) 
    return message


def send_email(info, imp):
    message = set_message(info, imp)
    server = 'smtp.gmail.com'
    sender = "r.soutelino@metocean.co.nz"
    smtp_client = smtplib.SMTP_SSL(host=server, port=465)
    smtp_client.login(sender, '@#Upwelling')
    receivers = ['r.soutelino@metocean.co.nz']
    smtp_client.sendmail(sender, receivers, message)


a = 10

try:
    while a < 20:
        if a == 10:
            a = 'Error'
        print a
        a += 1
        time.sleep(1)
except:
    send_email("FAILED", 'nsea')
    sys.exit(2)

send_email("Finished SUCCESSFULLY", 'nsea')