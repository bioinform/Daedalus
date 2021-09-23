import argparse
import pandas
import re

def umi_delim_extract(readName, regex):
    """
    Given a read name with UMI_R1 or UMI_R2 header labels from parse-umi, 
    extract the UMI information

    Parameters
    ----------
    readName : string
    The read header with umi labels, separated with ':'
    regex : string
    The UMI name: "UMI_R1" or "UMI_R2"
    
    Returns
    -------    
    the Position of UMI name in the read header, or None
    """
    stringVals = readName.split(":")    
    i = 0
    for val in stringVals:
        if re.search(regex, val):
            return(i)
        i += 1
    return None

def convert_quals(qualString):
    """
    Covert a string of quality integers '34,38,34,34...' to quality value characters 'CGCC' 

    Parameters
    ----------
    qualString : string
    read quality values integers, separated with a comma
    
    Returns
    -------    
    ASCII encoded quality value string
    """
    qualVals = []
    for qv in qualString.split(','):
        qualVals.append(chr(int(qv) + 33))
    return ''.join(qualVals)

def add_umi_columns(report_file, out_file):
    """
    Given a read report table, with read names labeled by `parse-umi`, 
    extract the UMI information and store the UMI string and quality values as columns.

    Parameters
    ----------
    report_file : path
    TSV report of reads with UMI labels
    out_file : path
    the name of the output tsv file
    
    Returns
    -------    
    None
    """    
    df = pandas.read_csv(report_file, sep="\t") 
    if not 'name' in df.columns:
        raise ValueError("read header column must be stored as 'name'")    
    ##extract UMI information
    checkName = df["name"][0]
    umi1Pos = umi_delim_extract(checkName, "UMI_R1")
    umi2Pos = umi_delim_extract(checkName, "UMI_R2")
    ##parse UMI data
    if umi1Pos is not None:
        df["umi1"] = df["name"].str.split(":", expand=True)[umi1Pos]
        df[['umi_r1_name', 'umi_r1', 'umi_r1_qual']] = df["umi1"].str.split("__", expand=True)
        df["umi_r1_qual"] = df["umi_r1_qual"].apply(convert_quals)
    if umi2Pos is not None:
        df["umi2"] = df["name"].str.split(":", expand=True)[umi2Pos]
        df[['umi_r2_name', 'umi_r2', 'umi_r2_qual']] = df["umi2"].str.split("__", expand=True)
        df["umi_r2_qual"] = df["umi_r2_qual"].apply(convert_quals)
    ##write output           
    df.to_csv(out_file, sep="\t", index=False)



def main():
    parser = argparse.ArgumentParser(
        description="Parse UMI information from read header and store UMI information in new columns")
    parser.add_argument("-r", "--reportFile",
                        required=True,
                        help="tsv read report with UMI labeled read headers")
    args = parser.parse_args()
    outFile = args.reportFile.replace(".tsv", "_umi.tsv")
    add_umi_columns(args.reportFile, outFile)

if __name__ == '__main__':
    main()
