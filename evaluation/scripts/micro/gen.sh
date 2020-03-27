echo 'cat ../input/1M.txt |' | tee 1.sh 10.sh 100.sh 1000.sh > /dev/null

for (( i = 0; i < 1; i++ )); do
  echo 'tr " " " " |' >> 1.sh
done

echo 'tr " " " " > out.txt' >> 1.sh

for (( i = 0; i < 10; i++ )); do
  echo 'tr " " " " |' >> 10.sh
done

echo 'tr " " " " > out.txt' >> 10.sh

for (( i = 0; i < 100; i++ )); do
  echo 'tr " " " " |' >> 100.sh
done

echo 'tr " " " " > out.txt' >> 100.sh

for (( i = 0; i < 1000; i++ )); do
  echo 'tr " " " " |' >> 1000.sh
done

echo 'tr " " " " > out.txt' >> 1000.sh
