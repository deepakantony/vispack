template <class T>
int write_raw(const char* name, const VISImage<T>& im)
{
    FILE* out_file;
    char filename[80];
    sprintf(filename, "%s", name);
    if ((out_file = fopen(filename, "w")) == NULL)
	{
	    printf("ERROR: output file open failed\n");
	    return(-1);
	}
    fwrite((const void*)(im.rep()->buffer()), sizeof(T),
	   im.height()*im.width(), out_file);
    return(0);
}

template <class T>
int read_raw(const char* name, VISImage<T>& im)
{
    FILE* out_file;
    char filename[80];
    if ((out_file = fopen(name, "w")) == NULL)
	{
	    printf("ERROR: output file open failed\n");
	    return(-1);
	}
    fread((void*)(im.repRef()->bufferRef()), sizeof(T),
	   im.height()*im.width(), 
	   out_file);
    return(0);
}
