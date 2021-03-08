import psycopg2
import datetime

try:
    connection = psycopg2.connect(user="xpduxdrvhsrgeh",
                                  password="74a8260bda6ab3f64fb2c7e6c8f02c48fb94a418d1e6c9d51356dcfffd93b2d3",
                                  host="ec2-18-211-97-89.compute-1.amazonaws.com",
                                  port="5432",
                                  database="d8vobvtjnnpf9j")
    cursor = connection.cursor()
    lines = tuple(open('experiment.txt', 'r'))
    postgres_insert_query = """ INSERT INTO Experiment 
		(experimentCode,username,contributeTime,assembleTime,calcstiffGlobalTime,solveGlobalTime,k,n,dim,refLevel,timeSelector,fechahora)
 		VALUES(%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)"""
    record_to_insert = (lines[0],lines[1],lines[2],lines[3],lines[4],lines[5],lines[6],lines[7],lines[8],lines[9],lines[10],'NOW()')
    cursor.execute(postgres_insert_query, record_to_insert)

    connection.commit()
    count = cursor.rowcount
    print(count, "Record inserted successfully into mobile table")

except (Exception, psycopg2.Error) as error:
    print("Failed to insert record into mobile table", error)

finally:
    # closing database connection.
    if connection:
        cursor.close()
        connection.close()
        print("PostgreSQL connection is closed")




