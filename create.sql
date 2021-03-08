SELECT * FROM experiment where k=1 AND n=2 AND dim=2 AND reflevel=4
SELECT * FROM experiment WHERE username LIKE '%ricardo-MAC%'

DELETE FROM experiment

CREATE TABLE experiment
(
    id integer NOT NULL GENERATED ALWAYS AS IDENTITY  (INCREMENT 1 ),
    username VARCHAR(50) NOT NULL,
    experimentcode VARCHAR(100)  NOT NULL,
    contributetime numeric(12,6) NOT NULL,
    assembletime numeric(12,6) NOT NULL,
    calcstiffglobaltime numeric(12,6) NOT NULL,
    solveglobaltime numeric(12,6) NOT NULL,
    k integer NOT NULL,
    n integer NOT NULL,
    dim integer NOT NULL,
    reflevel integer NOT NULL,
    timeselector VARCHAR(150) NOT NULL,
    fechahora timestamp NOT NULL
)
DROP TABLE experiment