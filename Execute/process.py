import psycopg2 as pg
import pandas.io.sql as psql
connection = pg.connect("host=ec2-18-211-97-89.compute-1.amazonaws.com dbname=d8vobvtjnnpf9j user=xpduxdrvhsrgeh password=74a8260bda6ab3f64fb2c7e6c8f02c48fb94a418d1e6c9d51356dcfffd93b2d3")
dataframe = psql.read_sql('SELECT solveglobaltime FROM Experiment WHERE reflevel=8', connection)
product_category = psql.read_sql_query('SELECT solveglobaltime FROM Experiment WHERE reflevel=7', connection)