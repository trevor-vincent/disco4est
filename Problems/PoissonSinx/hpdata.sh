sed -i 's/min_level = [0-9]/min_level = 0/g' options.input
./uniform_poisson_sinx_driver | grep D4EST_OUTPUT | cut -c 17- >> hp.dat
sed -i 's/min_level = [0-9]/min_level = 1/g' options.input
./uniform_poisson_sinx_driver | grep D4EST_OUTPUT | cut -c 17- >> hp.dat
sed -i 's/min_level = [0-9]/min_level = 2/g' options.input
./uniform_poisson_sinx_driver | grep D4EST_OUTPUT | cut -c 17- >> hp.dat
sed -i 's/min_level = [0-9]/min_level = 3/g' options.input
./uniform_poisson_sinx_driver | grep D4EST_OUTPUT | cut -c 17- >> hp.dat
sed -i 's/min_level = [0-9]/min_level = 4/g' options.input
./uniform_poisson_sinx_driver | grep D4EST_OUTPUT | cut -c 17- >> hp.dat
