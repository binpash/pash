import sqlite3
import os

connection = sqlite3.connect('database.db')

PASH_TOP=os.getenv('PASH_TOP')

with open(f'{PASH_TOP}/compiler/debug/schema.sql') as f:
    connection.executescript(f.read())

cur = connection.cursor()

cur.execute("INSERT INTO logs (worker, stderr) VALUES (?, ?)",
            ('worker1', 'Content for the first post')
            )

cur.execute("INSERT INTO logs (worker, stderr) VALUES (?, ?)",
            ('worker2', 'Content for the second post')
            )

connection.commit()
connection.close()