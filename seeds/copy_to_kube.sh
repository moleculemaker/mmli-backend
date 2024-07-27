kubectl exec -it mmli-backend-postgresql-0 -n alphasynthesis -- mkdir -p /opt/seeds/
for seed_file in ./*sql
do
    kubectl cp $seed_file mmli-backend-postgresql-0:/opt/seeds/ -n alphasynthesis
done
