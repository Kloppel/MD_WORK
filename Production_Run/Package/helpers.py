class helpers():
    def zip_data_with_labels(reduced, labels):
        rd_dcd = reduced[:, :end1]  
        rd_dcd2 = reduced[:, end1:(end1+add1)] 
        rd_namd = reduced[:,(end1+add1):]
        return zip([rd_dcd, rd_dcd2, rd_namd], labels)