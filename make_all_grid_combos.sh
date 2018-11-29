#!/bin/bash -e 

case *$HOSTNAME* in
        *ollie* )
                POOL_DIR_pism=/work/ollie/pgierz/pool_pism
                ;;
        *mlogin*)
                POOL_DIR_pism=/pf/a/a270077/pool_pism
                ;;
        *)
                echo "Unknown default PISM POOL location"
                exit
                ;;
esac

for dom in nhem greenland greenland_handbook laurentide; do
      for res in 5 10 20 40; do
              ./run_create_domain_grid.sh $dom $res ${POOL_DIR_pism}/grids
      done
done

for dom in antarctica; do
      for res in 4 8 16 32; do
              ./run_create_domain_grid.sh $dom $res ${POOL_DIR_pism}/grids
      done
done

echo "To create antarctica grids at 1 and 2 km resolution,"
echo "the ollie compute nodes run out of memory."
echo "You might need to ask for more memory or use a fat node."
echo
echo "Then run"
echo "   ./run_create_domain_grid.sh antarctica 1 ${POOL_DIR_pism}/grids"
echo "and"
echo "   ./run_create_domain_grid.sh antarctica 2 ${POOL_DIR_pism}/grids"
echo
