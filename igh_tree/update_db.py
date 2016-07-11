import sys
import time
import MySQLdb

infile = sys.argv[1]
table = sys.argv[2]
mysql_db_host = sys.argv[3]
mysql_db_user = sys.argv[4]
mysql_db_password = sys.argv[5]
mysql_db_name = sys.argv[6]

start_time = time.time()

# Open connection to database
db = MySQLdb.connect(host=mysql_db_host, # your host, usually localhost
                     user=mysql_db_user, # your username
                     passwd=mysql_db_password, # your password
                     db=mysql_db_name) # name of the data base
db_cursor = db.cursor()

def add_single_quotes_if_not_null(s):
    if s != "NULL":
        return "'" + s + "'"
    else:
        return s

# Update rows
with open(infile, 'rU') as f:
    
    col_names = f.readline().rstrip().split("\t")

    query_base = """
                 UPDATE {t} SET {x} WHERE uid = {u};
                 """

    for line in f:

        vals = line.rstrip().split("\t")
        uid = vals[0]
        key_val_pairs = ", ".join([n + " = " + add_single_quotes_if_not_null(v) for v, n in zip(vals[1:], col_names[1:])])

        subs = {"t": table, "x": key_val_pairs, "u": uid}

        query = query_base.format(**subs)

        try:
            db_cursor.execute(query)
        except MySQLdb.Error, e:
            print "Query failed"
            print query
            print e
            db.rollback()
            db.close()

print "Done!!"
print "update_db elapsed time (wall clock):", time.time() - start_time
