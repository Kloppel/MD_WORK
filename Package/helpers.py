class helpers():
    def zip_data_with_labels(reduced):
        rd_dcd = reduced[:, :98]  # first 98 frames
        rd_dcd2 = reduced[:, 98:(98+102)]  # next 102 frames
        rd_namd = reduced[:,(98+102):]  # last 100 frames
        return zip([rd_dcd, rd_dcd2, rd_namd], labels)