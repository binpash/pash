export FLASK_APP=$PASH_TOP/compiler/debug/app.py
export FLASK_ENV=development

# reset db
rm -f $PASH_TOP/compiler/debug/database.db
python3 $PASH_TOP/compiler/debug/init_db.py

# start server
flask run -p 5001 -h 0.0.0.0

# example way to run:
# docker exec nodemanager1 bash /opt/dish/pash/compiler/debug/server.sh
