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

# Insert rows
with open(infile, 'rU') as f:
    
    col_names = f.readline().rstrip().split("\t")
    col_names_formatted = ", ".join(col_names)
    field_dict = dict(zip(col_names, range(len(col_names))))

    select_query_base = """
                        SELECT EXISTS (
                        SELECT 1 FROM {t}
                        WHERE lineage_uid = '{lineage_uid}'
                        AND sequence_uid = '{sequence_uid}'
                        AND name = '{name}'
                        LIMIT 1 )
                        """

    update_query_base = """
                        UPDATE {t} SET value = '{value}'
                        WHERE lineage_uid = '{lineage_uid}'
                        AND sequence_uid = '{sequence_uid}'
                        AND name = '{name}';
                        """

    insert_query_base = """
                        INSERT INTO {t} (lineage_uid, sequence_uid, name, value)
                        VALUES ('{lineage_uid}', '{sequence_uid}', '{name}', '{value}');
                        """

    for line in f:

        vals = line.rstrip().split("\t")

        lineage_uid = vals[field_dict["lineage_uid"]]
        sequence_uid = vals[field_dict["sequence_uid"]]
        name = vals[field_dict["name"]]
        value = vals[field_dict["value"]]

        subs = {"t": table, "lineage_uid": lineage_uid, "sequence_uid": sequence_uid, "name": name, "value": value}

        select_query = select_query_base.format(**subs)
        update_query = update_query_base.format(**subs)
        insert_query = insert_query_base.format(**subs)

        # Test if row exists
        try:
            db_cursor.execute(select_query)
        except MySQLdb.Error, e:
            print "Query failed"
            print select_query
            print e
            db.rollback()
            db.close()
            
        if db_cursor.fetchone()[0]: # Exists

            try:
                db_cursor.execute(update_query)
            except MySQLdb.Error, e:
                print "Query failed"
                print update_query
                print e
                db.rollback()
                db.close()

        else:

            try:
                db_cursor.execute(insert_query)
            except MySQLdb.Error, e:
                print "Query failed"
                print insert_query
                print e
                db.rollback()
                db.close()

print "Done!!"
print "update_or_add_features elapsed time (wall clock):", time.time() - start_time
