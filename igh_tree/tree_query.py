import sys
import os
import time
import MySQLdb
import json

global verbose
verbose = False

def main(argv):

    ### Parse arguments
    infile_uids = argv[1]
    num_workers = int(argv[2])
    i_worker = int(argv[3].split(".")[-1])
    mysql_db_host = argv[4]
    mysql_db_user = argv[5]
    mysql_db_password = argv[6]
    mysql_db_name = argv[7]

    wdir = os.path.dirname(infile_uids)

    ### Get unique identifier of clustering
    with open(infile_uids, 'rU') as f:
        patient_uid = f.readline().rstrip()
        clustering_uid = f.readline().rstrip()

    ### Set output file name
    outfile_query = wdir + "/tree.query.temp." + str(i_worker)
    outfile_query_field_dict = wdir + "/tree.query.field_dict.json.temp." + str(i_worker)

    print "infile_uids", infile_uids
    print "outfile_query", outfile_query

    ### Query MySQL database for input
    start_time = time.time()

    field_dict = query(patient_uid, clustering_uid, num_workers, i_worker, outfile_query,
                       mysql_db_host, mysql_db_user, mysql_db_password, mysql_db_name)

    with open(outfile_query_field_dict, 'wb') as f:
        json.dump(field_dict, f)

    print "Query done!!"
    print "Elapsed time (wall clock):", time.time() - start_time

def query(patient_uid, clustering_uid, num_workers, i_worker, outfile,
          mysql_db_host, mysql_db_user, mysql_db_password, mysql_db_name):

    query = """
            SELECT l.uid, l.path
            FROM lineages l
            WHERE l.patient_uid = {p} AND l.clustering_uid = {c}
            AND {w} = {i}
            AND l.alignment = 1
            INTO OUTFILE '{f}'
            FIELDS TERMINATED BY '\\t'
            ENCLOSED BY ''
            LINES TERMINATED BY '\\n'
            """

    work_group = "l.uid_work_group_" + str(num_workers)

    subs = {"p": patient_uid, "c": clustering_uid, "w": work_group, "i": i_worker, "f": outfile}
    query = query.format(**subs)

    field_names = ["lineage_uid", "path"]
    field_dict = {k: i for i, k in enumerate(field_names)}

    try:
        os.remove(outfile)
    except OSError:
        pass

    print "Query is"
    print query

    start_time = time.time()

    db = MySQLdb.connect(host=mysql_db_host, # your host, usually localhost
                         user=mysql_db_user, # your username
                         passwd=mysql_db_password, # your password
                         db=mysql_db_name) # name of the data base
    db_cursor = db.cursor()

    print "Time to connect to db was"
    print time.time() - start_time

    start_time = time.time()

    try:
        db_cursor.execute(query)
    except MySQLdb.Error, e:
        print "Query failed"
        print query
        print e
        db.rollback()
        db.close()

    print "Query time was"
    print time.time() - start_time

    return field_dict

if __name__ == "__main__":
    main(sys.argv)
