# coding=utf-8

import matplotlib.pyplot as plt



# ATOM H1     HGA3    0.137 !    0.000
# ATOM H2     HGA3    0.123 !    0.000
# ATOM H3     HGA3    0.093 !    0.000
# ATOM S      SG311  -0.208 !   10.158
# ATOM C1     CG331  -0.25  !    2.500



patch1 = """

ATOM C2     CG321  -0.264 !   10.080
ATOM C3     CG2O1   0.714 !    6.966
ATOM O1     OG2D1  -0.544 !    1.805
ATOM N      NG2S1  -0.575 !    2.500
ATOM C4     CG2R61  0.29 !    0.000
ATOM C5     CG2R61 -0.225 !    0.000
ATOM C6     CG2R61 -0.179 !    1.899
ATOM C7     CG2R67 -0.028 !   21.091
ATOM C8     CG2R61 0.102 !   23.177
ATOM C9     CG2R61 -0.336 !    0.000
ATOM C10    CG2O3   0.581 !   20.266
ATOM O2     OG2D2  -0.607 !    0.025
ATOM O3     OG2D2  -0.69 !    0.025
ATOM C11    CG2R67 0.262 !   17.767
ATOM C12    CG2DC1 -0.02 !   25.770
ATOM C13    CG2DC1 -0.424 !   11.426
ATOM C14    CG2O5   0.662 !   26.350
ATOM C15    CG2DC1 -0.581 !  112.274
ATOM C16    CG2R61 -0.578 !    0.000
ATOM C17    CG2R61  0.661 !    0.000
ATOM C18    CG2R61 -0.423 !    0.000
ATOM C19    CG2R61 -0.029 !    0.000
ATOM C20    CG2R61  -0.222 !    1.595
ATOM C21    CG2R61  0.403 !   14.554
ATOM C22    CG2R61 -0.231 !   28.438
ATOM C23    CG2R61  0.409 !  109.816
ATOM O4     OG3R60 -0.306 !   79.888
ATOM O5     OG2D3  -0.665 !    2.506
ATOM O6     OG311  -0.665 !    0.000

ATOM H4     HGP1    0.291 !    0.000
ATOM H5     HGR61   0.156 !    0.000
ATOM H6     HGR61   0.089 !    0.000
ATOM H7     HGR61   0.145 !    0.000
ATOM H8     HGA4    0.097 !    1.425
ATOM H9     HGA4    0.113 !    2.673
ATOM H10    HGA4    0.15 !    3.483
ATOM H11    HGR61   0.149 !    0.000
ATOM H12    HGR61   0.113 !    0.000
ATOM H13    HGR61   0.102 !    0.000
ATOM H15    HGA2    0.136 !    0.000
ATOM H16    HGA2    0.112 !    0.000
"""

parse1 = patch1.split()
# print parse1
total = 0
k = 0
for i, entry in enumerate(parse1):
    if i == 3:
        total += float(entry)
        k = i + 6
    if i > 3 and i == k:
        total += float(entry)
        k += 6
print total

# nums_red = '0.36049968004226685 0.3392775356769562 0.3934885263442993 0.3867471516132355 0.35711050033569336 0.34814855456352234 0.38684025406837463 0.35138893127441406 0.36545541882514954 0.3636122941970825 0.3710690438747406 0.34680601954460144 0.30406489968299866 0.32542118430137634 0.3280116617679596 0.4100320637226105 0.36594247817993164 0.40875598788261414 0.4503396153450012 0.4695025384426117'
# reds = nums_red.split()
# red_rmsf = []
# red = open('/user/jdragelj/Desktop/red_rmsf_test.dat', 'w')
# for num in reds:
#     red.write(num+',\n')
#     red_rmsf.append(float(num))
# red.close()
#
# nums_oxi = '0.38689908385276794 0.397845983505249 0.44402366876602173 0.43040332198143005 0.4039624035358429 0.40566062927246094 0.4469844698905945 0.3964664340019226 0.38321301341056824 0.367086797952652 0.36816495656967163 0.36435842514038086 0.33901169896125793 0.3615415692329407 0.40429362654685974 0.4754070043563843 0.43557050824165344 0.419778972864151 0.4573356807231903 0.5094632506370544'
# oxis = nums_oxi.split()
# oxi_rmsf = []
# oxi = open('/user/jdragelj/Desktop/oxi_rmsf_test.dat', 'w')
# for num in oxis:
#     oxi.write(num+',\n')
#     oxi_rmsf.append(float(num))
# oxi.close()
#
# print len(oxis), len(reds)
#
# import numpy as np
# time_scale = np.arange(0.0,20.0,1.0)
# print time_scale
# print len(time_scale)
#
# fig = plt.figure()
# ax1 = fig.add_subplot(111)
# ax2 = ax1.twinx()
# ax1.plot(time_scale, red_rmsf, color='blue', lw=0.75)
# ax2.plot(time_scale, oxi_rmsf, color='red', lw=0.5)
# plt.show()