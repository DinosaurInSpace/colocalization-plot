from pandas import DataFrame as DF
import pickle

"""
The purpose of this script is to pull the errors from the hmdb_structure parser
log file "log.txt
"""

file_in = 'log.txt'

active_in = open(file_in, 'r')
count = 0
line_dict = {}
out_dict = {}
for line in active_in:
    count += 1
    line_dict[count] = line

    if 'ERROR' in line:
        previous = count - 1
        out_line = line_dict[previous]
        out_error = line_dict[count]
        out_tuple = (previous, count)
        out_list = [out_line, out_error]
        out_dict[out_tuple] = out_list

active_in.close()

error_df = DF.from_dict(out_dict, orient="index")
error_df.to_csv('error_log.csv', sep='\t')

outname = "error_df.pickle"
outfile = open(outname, "wb")
pickle.dump(error_df, outfile)
outfile.close()
