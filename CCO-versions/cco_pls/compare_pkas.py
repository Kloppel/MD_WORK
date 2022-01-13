# coding=utf-8


folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/Pr_xtra/prd/'
mds = ['md_Pr_prd-','md_Pr_prdh1','md_Pr_prdh2']
for md in mds:
    print md
    frames = [0,80,160,240,320,400,480]
    for frame in frames:

        file1 = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/Pr_xtra/prd/10ang_cavities_voda/%s/done/frame%i/results.dat' % (md, frame)
        file2 = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/Pr_xtra/prd/cavities_voda/%s/done/frame%i/results.dat' % (md, frame)

        f = open(file1, 'r')
        g = open(file2, 'r')

        fline = ''
        for line in f:
            if 'PRD-3_EHEM' in line:
                fline = line
        f.close()

        gline = ''
        for line in g:
            if 'PRD-3_EHEM' in line:
                gline = line
        g.close()

        fline = fline.split(':')
        gline = gline.split(':')

        print frame, float(fline[1].strip()) - float(gline[1].strip())
