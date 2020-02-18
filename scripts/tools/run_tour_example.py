import argparse
import csv
import os
import time
import pandas
parser = argparse.ArgumentParser(description= "Fakes running extended model")
parser.add_argument("-design","--design",dest="design",action="store",required=True)
parser.add_argument("-input","--input",dest="input",action="store",required=True)
parser.add_argument("-output","--output", dest="output", action="store", required=True)
args = parser.parse_args()

#args.design = os.path.abspath(args.design
args.output = os.path.abspath(args.output)
#args.input = os.path.abspath(args.input)

#time.sleep(10)

with open('bayesion_out.csv', 'wb') as file:
    writer = csv.writer(file)
    writer.writerow(["SN", "Name", "Contribution"])
    writer.writerow([1, "Linus Torvalds", "Linux Kernel"])
    writer.writerow([2, "Tim Berners-Lee", "World Wide Web"])
    writer.writerow([3, "Guido van Rossum", "Python Programming"])
df = pandas.read_csv('bayesion_out.csv')
df.to_csv(args.output, index=False)
