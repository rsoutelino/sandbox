#!/usr/bin/env python
from netCDF4 import Dataset

def nc3tonc4(filename4,filename3,clobber=False,nchunk=1000,quiet=False):
    """convert a netcdf 4 file (filename4) in NETCDF4_CLASSIC format
    to a netcdf 3 file (filename3) in NETCDF3_64BIT format."""
    ncfile4 = Dataset(filename4,'r')
    if ncfile4.file_format != 'NETCDF4_CLASSIC':
        raise IOError('input file must be in NETCDF4_CLASSIC format')
    ncfile3 = Dataset(filename3,'w',clobber=clobber,format='NETCDF3_64BIT')
    # create dimensions. Check for unlimited dim.
    unlimdimname = False
    unlimdim = None
    # create global attributes.
    if not quiet: print 'copying global attributes ..'
    for attname in ncfile4.ncattrs():
        setattr(ncfile3,attname,getattr(ncfile4,attname))
    if not quiet: print 'copying dimensions ..'
    for dimname,dim in ncfile4.dimensions.iteritems():
        if dim.isunlimited():
            unlimdimname = dimname
            unlimdim = dim
            ncfile3.createDimension(dimname,None)
        else:
            ncfile3.createDimension(dimname,len(dim))
    # create variables.
    for varname,ncvar in ncfile4.variables.iteritems():
        if not quiet: print 'copying variable',varname
        # is there an unlimited dimension?
        if unlimdimname and unlimdimname in ncvar.dimensions:
            hasunlimdim = True
        else:
            hasunlimdim = False
        var = ncfile3.createVariable(varname,ncvar.dtype,ncvar.dimensions)
        # fill variable attributes.
        for attname in ncvar.ncattrs():
            setattr(var,attname,getattr(ncvar,attname))
        # fill variables with data.
        if hasunlimdim: # has an unlim dim, loop over unlim dim index.
            # range to copy
            if nchunk:
                start = 0; stop = len(unlimdim); step = nchunk
                if step < 1: step = 1
                for n in range(start, stop, step):
                    nmax = n+nchunk
                    if nmax > len(unlimdim): nmax=len(unlimdim)
                    var[n:nmax] = ncvar[n:nmax]
            else:
		var[0:len(unlimdim)] = ncvar[:]
        else: # no unlim dim or 1-d variable, just copy all data at once.
            var[:] = ncvar[:]
        ncfile3.sync() # flush data to disk
    # close files.
    ncfile3.close()
    ncfile4.close()

if __name__ == '__main__':

    import sys, getopt, os

    usage = """
 Convert a netCDF 4 file (in NETCDF4_CLASSIC format) to netCDF 3 format.

 usage: %s [-h] [-o] [--chunk] netcdf4filename netcdf3filename
 -h -- Print usage message.
 -o -- Overwite destination file (default is to raise an error if output file already exists).
 --quiet=(0|1)  -- if 1, don't print diagnostic information.
 --chunk=(integer) -- number of records along unlimited dimension to 
     write at once.  Default 1000.  Ignored if there is no unlimited 
     dimension.  chunk=0 means write all the data at once.
\n""" % os.path.basename(sys.argv[0])

    try:
        opts, pargs = getopt.getopt(sys.argv[1:], 'ho',
                                    ['chunk=','quiet='])
    except:
        (type, value, traceback) = sys.exc_info()
        print "Error parsing the options. The error was:", value
        sys.stderr.write(usage)
        sys.exit(0)

    # default options
    quiet = 0
    chunk = 1000
    overwritefile = 0

    # Get the options
    for option in opts:
        if option[0] == '-h':
            sys.stderr.write(usage)
            sys.exit(0)
        elif option[0] == '-o':
            overwritefile = 1
        elif option[0] == '--quiet':
            quiet = int(option[1])
        elif option[0] == '--chunk':
            chunk = int(option[1])
        else:
            print option[0], ": Unrecognized option"
            sys.stderr.write(usage)
            sys.exit(0)

    # if we pass a number of files different from 2, abort
    if len(pargs) <> 2:
        print "You need to pass both source and destination!."
        sys.stderr.write(usage)
        sys.exit(0)

    # Catch the files passed as the last arguments
    filename4 = pargs[0]
    filename3 = pargs[1]

    # copy the data from filename4 to filename3.
    nc3tonc4(filename4,filename3,clobber=overwritefile,quiet=quiet)

