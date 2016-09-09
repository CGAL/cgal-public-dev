# cgal-public-dev
This is the Git repository used to store some on-going work by CGAL developers. The repository that hosts the `master` branch of CGAL is [`CGAL/cgal`](http://github.com/CGAL/cgal).

# Start working
To start working with this repository, you must also set the a remote repository that contains the `master` branch of CGAL:

    git remote add cgal git@github.com:CGAL/cgal.git

Then you can create your own working branch. Say, your name is *rineau* and you will work on a new read/write function for the `Triangulation` package:

    git checkout -b gsoc2015-Triangulation-add_input_output-rineau cgal/master
