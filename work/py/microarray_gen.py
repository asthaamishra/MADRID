#!/usr/bin/python3

import argparse
import re
import os
import sys
import pandas as pd
import numpy as np
from scipy import stats
from sqlalchemy import Column, ForeignKey, Integer, String, Float
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship, sessionmaker, load_only
from sqlalchemy import create_engine
from project import configs
from GSEpipelineFast import *


# create declarative_base instance
Base = declarative_base()

# creates a create_engine instance at the bottom of the file
engine = create_engine("sqlite:///microarray.db")


# Add SQL table information
class Sample(Base):
    __tablename__ = "sample"

    id = Column(Integer, primary_key=True)
    ENTREZ_GENE_ID = Column(String(250), nullable=False)
    VALUE = Column(Float, nullable=False)
    P_VALUE = Column(Float)
    ABS_CALL = Column(String(2))
    Sample = Column(String(64))


class IDMapps(Base):
    __tablename__ = "id_entrez_map"

    idx = Column(Integer, primary_key=True)
    ID = Column(String(250))
    ENTREZ = Column(String(250))


class GSEinfo(Base):
    __tablename__ = "gseinfo"

    Sample = Column(String(64), primary_key=True)
    GSE = Column(String(64))
    Platform = Column(String(64))


Base.metadata.create_all(engine)

DBSession = sessionmaker(bind=engine)
session = DBSession()
# --- Finish initializing database ---


def lookupMicroarrayDB(gseXXX):
    """
    check if gse already in database
    :param gseXXX:
    :return:
    """
    df = pd.read_sql(
        session.query(GSEinfo)
        .filter_by(GSE=gseXXX.gsename)
        .options(load_only("Sample", "GSE"))
        .statement,
        session.bind,
    )
    gsm_list = gseXXX.gsm_platform.keys()
    gsm_db = df.Sample.tolist()
    if set(gsm_list).issubset(set(gsm_db)):
        return True
    else:
        return False


def updateMicroarrayDB(gseXXX):
    """
    Update GSE info to microarray.db
    :param gseXXX: gse object
    :return:
    """
    df_clean = gseXXX.get_entrez_table_pipeline()
    if df_clean.empty:
        return None

    df_clean.sort_index(inplace=True)
    df_clean.reset_index(inplace=True)

    # write to database table sample
    df_samples = pd.DataFrame(
        [], columns=["ENTREZ_GENE_ID", "VALUE", "P_VALUE", "ABS_CALL", "Sample"]
    )
    df_samples.index.name = "id"

    cols_clean = list(df_clean)
    for key, val in gseXXX.gsm_platform.items():
        col_val = "{}.cel.gz".format(key.lower())
        col_abs = "{}.cel.gz.1".format(key.lower())
        col_p = "{}.cel.gz.2".format(key.lower())
        if col_val not in cols_clean:
            continue
        df_s = pd.DataFrame([])
        df_s["ENTREZ_GENE_ID"] = df_clean["ENTREZ_GENE_ID"]
        df_s["VALUE"] = df_clean[col_val]
        df_s["ABS_CALL"] = df_clean[col_abs]
        df_s["P_VALUE"] = df_clean[col_p]
        df_s["Sample"] = key.upper()
        df_s.index.name = "id"
        df_samples = pd.concat(
            [df_samples, df_s.dropna(how="any")], ignore_index=True, sort=True
        )

    df_samples.index.name = "id"
    df_samples.to_sql(con=engine, name="sample", if_exists="append", index=False)

    # write to database table GSEinfo
    df_gseinfo = pd.DataFrame([], columns=["Sample", "GSE", "Platform"])
    df_gseinfo["Sample"] = pd.Series(list(gseXXX.gsm_platform.keys()))
    df_gseinfo["Platform"] = pd.Series(list(gseXXX.gsm_platform.values()))
    df_gseinfo["GSE"] = gseXXX.gsename

    df_gseinfo.to_sql(con=engine, name="gseinfo", if_exists="append", index=False)


# function to complete the inquery of a sheet
def queryTest(df, expression_proportion, top_proportion):
    sr = df["GSE ID"]
    gse_ids = sr[sr.str.match("GSE")].unique()
    sr = df["Samples"].dropna()
    gsm_ids = sr.unique()
    print("---\nStart Collecting Data for:")
    print(gse_ids)
    print(gsm_ids)
    print("---\n")
    # fetch data of each gse if it is not in the database, update database
    for gse_id in gse_ids:
        querytable = df[df["GSE ID"] == gse_id]
        gseXXX = GSEproject(gse_id, querytable, configs.rootdir)
        if lookupMicroarrayDB(gseXXX):
            print("{} already in database, skip over.".format(gseXXX.gsename))
            continue
        updateMicroarrayDB(gseXXX)

    df_results = fetchLogicalTable(gsm_ids)
    df_output = mergeLogicalTable(df_results)

    df_output = df_output.apply(pd.to_numeric)
    posratio = df_output.sum(axis=1, skipna=True) / df_output.count(axis=1)

    df_output["Pos"] = posratio
    df_output["expressed"] = np.where(posratio >= expression_proportion, 1, 0)
    df_output["top"] = np.where(posratio >= top_proportion, 1, 0)

    return df_output


def fetchLogicalTable(gsm_ids):
    """
    Fetch the Logical Table Based on ABS_CALL of Samples
    :param gsm_ids: list of sample names to fetch
    :return: pandas dataframe
    """
    df_results = pd.DataFrame([], columns=["ENTREZ_GENE_ID"])
    df_results.set_index("ENTREZ_GENE_ID", inplace=True)
    
    for gsm in gsm_ids:
        df = pd.read_sql(
            session.query(Sample)
            .filter_by(Sample=gsm)
            .options(load_only("ABS_CALL", "ENTREZ_GENE_ID"))
            .statement,
            session.bind,
            index_col="ENTREZ_GENE_ID",
        )
        df.drop("id", axis=1, inplace=True)
        df.rename(columns={"ABS_CALL": gsm}, inplace=True)
        df.loc[df[gsm] == "1", gsm] = "A"
        df.loc[df[gsm] == "3", gsm] = "P"
        df.loc[df[gsm] == "2", gsm] = "M"

        df.loc[df[gsm] == "A", gsm] = 0
        df.loc[df[gsm] == "P", gsm] = 1
        df.loc[df[gsm] == "M", gsm] = 1
        df.sort_values(by=["ENTREZ_GENE_ID", gsm], inplace=True)
        df = df[~df.index.duplicated(keep="last")]

        df_results = pd.concat([df_results, df], axis=1, sort=False)

    # Need to set index name after merge
    df_results.index.name = "ENTREZ_GENE_ID"
    
    return df_results


# Merge Output
def mergeLogicalTable(df_results):
    """
    Merge the Rows of Logical Table belongs to the same ENTREZ_GENE_ID
    :param df_results:
    :return: pandas dataframe of merged table
    """
    # step 1: get all plural ENTREZ_GENE_IDs in the input table, extract unique IDs
    df_results.reset_index(drop=False, inplace=True)
    df_results["ENTREZ_GENE_ID"] = df_results["ENTREZ_GENE_ID"].astype(str)
    df_results["ENTREZ_GENE_ID"] = df_results["ENTREZ_GENE_ID"].str.replace(
        " /// ", "//"
    )
    id_list = []
    df_results.dropna(axis=0, subset=["ENTREZ_GENE_ID"], inplace=True)
    entrez_single_id_list = df_results[
        ~df_results["ENTREZ_GENE_ID"].str.contains("//")
    ]["ENTREZ_GENE_ID"].tolist()
    entrez_id_list = df_results[df_results["ENTREZ_GENE_ID"].str.contains("//")][
        "ENTREZ_GENE_ID"
    ].tolist()
    for entrez_id in entrez_id_list:
        entrez_ids = entrez_id.split("//")
        id_list.extend(entrez_ids)
        df_dups = pd.DataFrame(
            [], columns=list(df_results), index=list(range(len(entrez_ids)))
        )
        dup_rows = pd.DataFrame([])
        for eid in entrez_ids:
            rows = df_results.loc[df_results["ENTREZ_GENE_ID"] == entrez_id].copy()
            rows["ENTREZ_GENE_ID"] = eid
            dup_rows = pd.concat([dup_rows, rows], axis=0)
        df_results = pd.concat(
            [df_results, pd.DataFrame(dup_rows)], axis=0, ignore_index=True
        )

        df_results.drop(
            df_results[df_results["ENTREZ_GENE_ID"] == entrez_id].index, inplace=True
        )

    common_elements = list(set(entrez_single_id_list).intersection(set(id_list)))
    dups = [x for x in id_list if id_list.count(x) > 1]

    full_entre_id_sets = []
    cnt = 0
    entrez_dups_list = []
    idx_list = list(range(len(entrez_id_list)))
    
    for idx1 in range(len(entrez_id_list)):
        if idx1 not in idx_list:
            continue
            
        set1 = set(entrez_id_list[idx1].split("//"))
        idx_list.remove(idx1)
        toremove = []
        
        for idx2 in idx_list:
            set2 = set(entrez_id_list[idx2].split("//"))
            intersect = set1.intersection(set2)       
            if bool(intersect):
                set1 = set1.union(set2)
                toremove.append(idx2)
                
        for idx3 in toremove:
            idx_list.remove(idx3)
            
        sortlist = list(set1)
        sortlist.sort(key=int)
        new_entrez_id = " /// ".join(sortlist)
        full_entre_id_sets.append(new_entrez_id)
        
    full_entre_id_sets = list(set(full_entre_id_sets))

    for full_entrez_id in full_entre_id_sets:
        singles = full_entrez_id.split(" /// ")
        entrez_dups_list.append(singles)
        cnt += 1

    entrez_dups_dict = dict(zip(full_entre_id_sets, entrez_dups_list))

    for merged_entrez_id, entrez_dups_list in entrez_dups_dict.items():
        df_results["ENTREZ_GENE_ID"].replace(to_replace=entrez_dups_list, value=merged_entrez_id, inplace=True)

    df_results.set_index("ENTREZ_GENE_ID", inplace=True)
    df_output = df_results.fillna(-1).groupby(level=0).max()
    df_output.replace(-1, np.nan, inplace=True)

    # TODO: Test if this is working properly
    """
    There seems to be an error when running Step 2.1 in the pipeline.ipynb file
    The commented-out return statement tries to return the df_output dataframe values as integers, but NaN values exist
        Because of this, it is unable to do so.
    If we change this to simply output the database, the line "np.where(posratio >= top_proportion . . ." (line ~162)
        Fails because it is comparing floats and strings
    
    I am unsure what to do in this situation
    """
    
    #return df_output.astype(int)
    return df_output


def load_microarray_tests(filename, context_name):
    def load_empty_dict():
        savepath = os.path.join(
            configs.rootdir,
            "data",
            "data_matrices",
            "dummy",
            "dummy_microarray_data.csv",
        )
        dat = pd.read_csv(savepath, index_col="ENTREZ_GENE_ID")
        
        return "dummy", dat

    if (not filename or filename == "None"):  # if not using microarray use empty dummy matrix
    
        return load_empty_dict()

    inquiry_full_path = os.path.join(configs.rootdir, "data", "config_sheets", filename)
    if not os.path.isfile(inquiry_full_path): 
        print("Error: file not found {}".format(inquiry_full_path))
        sys.exit()

    filename = "Microarray_{}.csv".format(context_name)
    fullsavepath = os.path.join(
        configs.rootdir, "data", "results", context_name, "microarray", filename
    )
    if os.path.isfile(fullsavepath):
        data = pd.read_csv(fullsavepath, index_col="ENTREZ_GENE_ID")
        print("Read from {}".format(fullsavepath))
        
        return context_name, data
        
    else:
        print(
            f"Microarray gene expression file for {context_name} was not found at {fullsavepath}. This may be "
            f"intentional. Contexts where microarray data can be found in /work/data/results/{context_name}/ will "
            f"still be used for other contexts if found."
        )
        
        return load_empty_dict()


def main(argv):
    inputfile = "microarray_data_inputs.xlsx"

    parser = argparse.ArgumentParser(
        prog="microarray_gen.py",
        description="This file is for processing microarray data",
        epilog="For additional help, please post questions/issues in the MADRID GitHub repo at "
        "https://github.com/HelikarLab/MADRID or email babessell@gmail.com",
    )
    parser.add_argument(
        "-c",
        "--config-file",
        type=str,
        required=True,
        dest="config_file",
        help="The microarray configuration file name",
    )
    parser.add_argument(
        "-e",
        "--expression-proportion",
        type=float,
        required=True,
        dest="expression_proportion",
        help="Number of sources with active gene for it to be considered active even if it is not a high confidence-gene",
    )
    parser.add_argument(
        "-t",
        "--top-proportion",
        type=str,
        required=True,
        dest="top_proportion",
        help="Genes can be considered high confidence if they are expressed in a high proportion of samples. "
             + "High confidence genes will be considered expressed regardless of agreement with other data sources",
    )
    args = parser.parse_args()
    inputfile = args.config_file
    expression_proportion = args.expression_proportion
    top_proportion = args.top_proportion

    print("Input file is ", inputfile)
    print("Expression Proportion for Gene Expression is ", expression_proportion)
    print("Top proportion for high-confidence genes is ", top_proportion)

    inqueryFullPath = os.path.join(configs.rootdir, "data", "config_sheets", inputfile)
    xl = pd.ExcelFile(inqueryFullPath)
    sheet_names = xl.sheet_names
    inqueries = pd.read_excel(inqueryFullPath, sheet_name=sheet_names, header=0)

    for context_name in sheet_names:
        inqueries[context_name].fillna(method="ffill", inplace=True)
        df = inqueries[context_name].loc[:, ["GSE ID", "Samples", "GPL ID", "Instrument"]]
        df_output = queryTest(df, expression_proportion, top_proportion)
        filename = "Microarray_{}.csv".format(context_name)
        fullsavepath = os.path.join(
            configs.rootdir, "data", "results", context_name, filename
        )
        os.makedirs(os.path.dirname(fullsavepath), exist_ok=True)
        df_output.to_csv(fullsavepath)
        print("Save to {}".format(fullsavepath))


if __name__ == "__main__":
    main(sys.argv[1:])