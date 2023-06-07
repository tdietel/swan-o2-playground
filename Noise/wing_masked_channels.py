#!/usr/bin/env python3

import cx_Oracle




connection = cx_Oracle.connect(user='trdfero', password='trd1234', dsn='alice_fero')

sql = """
SELECT sm_slot,sm_stack,sm_layer,roc_type,roc_serial,roc_pos,rob_pos,adcmask 
FROM cferr_adc_current 
	JOIN rob_loc USING (rob_id) 
	JOIN roc_loc USING(roc_type,roc_serial) 
	JOIN sm_loc USING(sm_id)
ORDER BY sm_slot,sm_stack,sm_layer,roc_pos,rob_pos"""

print("Get all rows via iterator")
cursor = connection.cursor()

with open("masked_channels.csv", "w") as f:
    f.write("sm_slot,sm_stack,sm_layer,roc_type,roc_serial,roc_pos,rob_pos,adcmask\n")

    for result in cursor.execute(sql):
        f.write(",".join([str(x) for x in result])+"\n")
        print(result)
    print()

#print("Query one row at a time")
#cursor.execute(sql)
#row = cursor.fetchone()
#print(row)
#row = cursor.fetchone()
#print(row)
#print()

#print("Fetch many rows")
#cursor.execute(sql)
#res = cursor.fetchmany(numRows=3)
#print(res)

