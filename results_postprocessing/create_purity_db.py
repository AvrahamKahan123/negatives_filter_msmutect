import os
import sqlite3
from dataclasses import dataclass
from typing import List


@dataclass
class Locus:
    chromosome: str
    start: str
    stop: str
    pattern: str
    purity: float


class DB_connection:
    def __init__(self, sqlite_fp: str):
        self.connection = sqlite3.connect(sqlite_fp)
        self.cursor = self.connection.cursor()
        self.table_name = "purity"
        self.encoded_location_column_name = "encoded_location"
        self.purity_column_name = "purity"

    @staticmethod
    def db_path():
        return "C:/Users/avrah/MaruvkaLab/post_processing_code_for_australians/purity.db"

    def create_purity_table(self):
        creation_statement = f"""
            CREATE TABLE IF NOT EXISTS {self.table_name} (
                {self.encoded_location_column_name} TEXT PRIMARY KEY,
                {self.purity_column_name} FLOAT NOT NULL
            )
            """
        self.cursor.execute(creation_statement)

    def add_rows(self, loci: List[Locus]):
        encoded_loci = [self.encode_locus(locus) for locus in loci]
        rows_to_add = []
        for encoded_locus, locus in zip(encoded_loci, loci):
            rows_to_add.append((encoded_locus, locus.purity))
        self.cursor.executemany(f"INSERT INTO {self.table_name} ({self.encoded_location_column_name}, {self.purity_column_name}) VALUES (?, ?)", rows_to_add)

    def encode_locus(self, locus: Locus):
        return f"{locus.chromosome}_{locus.start}_{locus.stop}_{locus.pattern}"

    def query_locus_purity(self, locus: Locus) -> float:
        encoded_locus = self.encode_locus(locus)
        self.cursor.execute(f"SELECT * FROM purity WHERE {self.encoded_location_column_name} = ?", (encoded_locus,))
        rows = self.cursor.fetchall()
        try:
            return float(rows[0][1])
        except IndexError:
            raise IndexError(f"Could not find Locus: {locus}")

    def close(self):
        self.connection.commit()
        self.connection.close()

    def __del__(self):
        self.close()

def connect_to_purity_db():
    return DB_connection(DB_connection.db_path())


def main(loci_fp="C:/Users/avrah/MaruvkaLab/GRCh38.d1.vd1_1to15_repetitive_loci_sorted_fixed.phobos"):
    # print(os.path.exists(DB_connection.db_path()))
    # with open(loci_fp, 'r') as loci_f:
    #     rows = []
    #     # i=0
    #     for line in loci_f:
    #         # i+=1
    #         split_line = line.split("\t")
    #         rows.append(Locus(split_line[0], split_line[3], split_line[4], split_line[12], round(float(split_line[2]))))
    #         print(rows[-1])
    #         exit()
            # if i==10:
            #     break
    # db = connect_to_purity_db()
    # db.create_purity_table()
    # db.add_rows(rows)
    db = connect_to_purity_db()

    print(db.query_locus_purity(Locus("chr1", str(10001), str(10468), "AACCCT", None)))



if __name__ == '__main__':
    main()
