#!/usr/bin/env python

import argparse
import pandas as pd
import os.path
import smtplib

from os import path

def parse_args():
    parser = argparse.ArgumentParser(description="""Send email when job completes""")
    parser.add_argument('--email', type=str, help='Email address to send pipeline run summary to')
    parser.add_argument('--summary', type=str, help='Path to pipeline output summary')
    parser.add_argument('--name', type=str, help='Pipeline run name')
    args = parser.parse_args()
    return(args)

def send_email(message, email):
    sender = 'erin.kesel@contractors.roche.com'
    receiver = email
    
    try:
        smtpObj = smtplib.SMTP('localhost')
        smtpObj.sendmail(sender, receiver, message)
    except SMTPException:
        print ("Error: unable to send email")

def success_message(samplesum, jobname):
    col_list = ["status"]
    summary = pd.read_csv(samplesum, usecols=col_list)
    message = "Your Daedalus pipeline run %s has finished successfully.\n\nRun summary:\n" % jobname
    message += "You submitted %i samples\n" % summary.shape[0]
    for i, v in summary["status"].value_counts().items():
        message += '%i samples %s\n' % (v,i)
    return(message)

def fail_message():
    message = "Your Daedalus pipeline run has failed. Please contact the developer.\n"

def main(args):
    subject = "Daedalus pipeline run - %s" %args.name
    if(path.exists(args.summary)): 
        body = success_message(args.summary, args.name)
    else:
        body = fail_message()
    message = 'Subject: {}\n\n{}'.format(subject, body)
    send_email(message, args.email)

if __name__ == '__main__':
    args = parse_args()
    main(args)


