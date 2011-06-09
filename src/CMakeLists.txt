# Recurse into the "Hello" and "Demo" subdirectories.  This does not actually
# cause another cmake executable to run.  The same process will walk through
# the project's entire directory structure.

#this is a comment line

ADD_DEFINITIONS( -g )
SUBDIRS(util image mat vol scan voxmodel diffuse param)

#ADD_DEFINITIONS( -g )