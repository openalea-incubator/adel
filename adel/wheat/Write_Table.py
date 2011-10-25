def Write_Table(Table,filename,first=True):
    '''
    returns a string representation of the content of the table
    '''
    
    header = [' '.join(Table.keys())]
    rows = zip(*[Table[k] for k in Table.keys()])
    lines = [' '.join(map(str,r)) for r in rows]
    mode = "a"
    if first:
        mode = "w"
    fout = open(filename,mode)
    if first :
        fout.write('\n'.join(header + lines))
    else :
        fout.write('\n')
        fout.write('\n'.join(lines))
    fout.close()
    
    return filename,
