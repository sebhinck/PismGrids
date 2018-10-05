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

echo "NOTE: antarctica resolutions are slightly different, these are not yet made by default."
echo "NOTE: if you are using antarctica and need a griddes file at a resolution that isnt"
echo "NOTE: available, please contact the current pism_pool maintainer (Paul Gierz)"

