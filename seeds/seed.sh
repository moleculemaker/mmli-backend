for seed_file in ./*.sql
do
  psql -U postgres -d mmli -f $seed_file
done