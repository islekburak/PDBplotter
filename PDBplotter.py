import argparse, textwrap
from collections import OrderedDict
import pandas as pd


parser=argparse.ArgumentParser(prog="PDBplotter.py",
	usage="python3 %(prog)s [options] <path_of_file>",
	description=textwrap.dedent("""\
		authors:	islekbro
		contact:	https://github.com/islekburak"""),
	add_help=True,
	formatter_class=argparse.RawTextHelpFormatter, 
	epilog=textwrap.dedent("""\

#To use this program, you should first create a new file from PDB file using VMD. To do this, please follow the steps below:
===========================================================================================================================================

# open your PDB molecule in VMD
# open Tk Console (in extensions tab) and write:

#	<atomselect keywords> 	(to see the list of keywords)
#	<set all [atomselect top "protein"]>
#	<atomselect0 get structure>
#	<atomselect0 get resid>
#	<atomselect0 get resname>

# copy & paste all data into a .txt file line by line (via disable word wrap function in Sublime Text). Save this file as whatever.txt
# Now you will use this whatever.txt file to parse data and create a plot.
===========================================================================================================================================
"""))

parser.add_argument("-i","--input", metavar="<.TXT FILE>", required=True, help="Write the path of your input -txt- file (i.e. /home/Desktop/input.txt)")
parser.add_argument("-o","--outdir", metavar="<OUTDIR>", required=True, help="Write only the path where you want to put outs (i.e. /home/Desktop/)")
args=parser.parse_args()

#getting directory
directory = args.outdir

#getting input path
inp_file=args.input
inpname=inp_file.split("/")[-1].split(".")[0]

a=pd.read_csv(inp_file,sep=" ",header=None)
a=a.transpose()

headers=["structure","resid","resname"]
a=pd.DataFrame(a.values, columns =headers)

a["key"]= (a["structure"] != a['structure'].shift(1)).astype(int).cumsum()


b=a.groupby(['key','structure'])['resid'].apply(','.join)
c=a.groupby(['key','structure'])['resname'].apply(','.join)

df1=pd.DataFrame(b)
df1=df1.reset_index()
df2=pd.DataFrame(c)
df2=df2.reset_index()


df1["resname"]=df2["resname"]
df1["resid"]=df1["resid"].astype(str)

listo=[]
for i in range(len(df1)):
	val=df1["resid"][i].split(",")[0]+","+df1["resid"][i].split(",")[-1]
	listo.append(val)

df1["range"]=listo
df1.drop(columns =["resid"], inplace = True)

new=df1["range"].str.split(",",expand=True)
df1["resid_start"]= new[0]
df1["resid_end"]= new[1]
df1.drop(columns =["range"], inplace = True)
df1.drop(columns =["key"], inplace = True)

reso=[]
for i in range(len(df1)):
	burak=",".join(OrderedDict.fromkeys(df1["resname"][i].split(',')))
	reso.append(burak)

df1["residues"]=reso
df1.drop(columns =["resname"], inplace = True)
df1["protein"]=inpname

cols=list(df1.columns)
cols=[cols[-1]]+cols[:-1]
df=df1[cols]

mergedpath = directory+ "/"+inpname+"_PDBparsed.csv"
df.to_csv(mergedpath, index=None)