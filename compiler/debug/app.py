import sqlite3
from flask import Flask, render_template, request, url_for, flash, redirect, jsonify
from werkzeug.exceptions import abort

app = Flask(__name__)
app.config['SECRET_KEY'] = 'my very secret key that is not used'
worker_db = {}

def get_db_connection():
    conn = sqlite3.connect('database.db')
    conn.row_factory = sqlite3.Row
    return conn

def get_logs(worker):
    conn = get_db_connection()
    worker_logs = conn.execute('SELECT * FROM logs WHERE worker = ?',
                        (worker,)).fetchall()
    conn.close()
    if worker_logs is None:
        abort(404)
    return worker_logs

@app.route('/')
def index():
    conn = get_db_connection()
    workers = conn.execute('SELECT DISTINCT worker FROM logs').fetchall()
    conn.close()
    return render_template('index.html', workers=workers)

@app.route('/fail')
def fail():
    conn = get_db_connection()
    logs = conn.execute('SELECT * FROM logs WHERE returncode != 0').fetchall()
    conn.close()
    
    if logs is None:
        abort(404)
    return render_template('fail.html', logs=logs)

@app.route('/worker/<string:worker>')
def worker_log(worker):
    worker_logs = get_logs(worker)
    return render_template('wlog.html', logs=worker_logs)

@app.route('/debug', methods=['POST'])
def debug():
    for key, val in request.json:
        print(key, val.stderr.replace("\n", "<br />"))

@app.route('/putlog', methods=['POST'])
def put_log():
    worker = request.json['name']
    stderr = request.json['stderr']
    returncode = request.json['returncode']
    shellscript = request.json['shellscript']
    
    if not worker:
        abort(401)
    
    stderr = stderr.replace("\n", "<br />")
    shellscript = shellscript.replace("\n", "<br />")

    print(worker)
    print(stderr)
    # print(shellscript)
    print()

    conn = get_db_connection()
    conn.execute('INSERT INTO logs (worker, returncode, stderr, shellscript) VALUES (?, ?, ?, ?)',
                         (worker, returncode, stderr, shellscript))
    
    conn.commit()
    conn.close()
    return jsonify({'code': 201})

@app.route('/clearall', methods=['POST'])
def clearall():    
    conn = get_db_connection()
    c = conn.cursor()
    c.execute('DELETE FROM logs')
    
    conn.commit()
    conn.close()
    return jsonify({'code': 201})

class Worker:
    def __init__(self, id):
        self.id = id
        self.log = ""

    def clear_log(self):
        self.log = ""