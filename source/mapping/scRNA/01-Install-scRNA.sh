#!/bin/bash

cd /projects/compbio/users/xran2/tool/scRNA_map

#-------------------------------------------#
#                                           #
#           Part 1: Download Files          #
#                                           #
#-------------------------------------------#

# Download Cell Ranger 9.0.1 (Feb 6, 2025)
# Link: https://www.10xgenomics.com/support/software/cell-ranger/downloads
# Cell Ranger 9.0.1 (Feb 6, 2025)
# File size: 807 MB
# md5sum: md5sum: 2efec98bff01f7a59edaf43724fae13f

curl -o cellranger-9.0.1.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-9.0.1.tar.gz?Expires=1738944703&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=mvHcCGkjl~4RtiOBsZK-pQ62Tkgj9HhYnynNGjQvQ0IixB9JxUihb4P0j3Aj4UPj3ji-ClBzR1i07NUwyd~6M1kVfFuF3TqPtqNMn8jh2RhEJhy3BrEyAcmp~bjKLcGGYTkLLMl-7j558eTphlgExCM~0679xn~CUnC4HNGwNG7DhsDrRGbdI~p3kXYtDfwWP1r1BZM9s1oYYxXYh9Yv5V8Ty9J6dmoW5xVfnlp9atZKjvRHqNj2Zmp8sHFVkjRNQ0~Yo4rX20FwJmFU8OObUu8hQ12BkcEIFIJRzHdBRoQnTjK9kSmtdQp7yR71VfrCYrkmcGsTzftZ27ZamfdbeQ__"

# Check the md5sum
md5sum cellranger-9.0.1.tar.gz

# Download Human reference data
# Link: https://www.10xgenomics.com/support/software/cell-ranger/downloads
# Human reference (GRCh38) - 2024-A
# File size: 11 GB
# md5sum: a7b5b7ceefe10e435719edc1a8b8b2fa

curl -O "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz"

# Check the md5sum
md5sum refdata-gex-GRCh38-2024-A.tar.gz


#-------------------------------------------#
#                                           #
#           Part 2: Installation            #
#                                           #
#-------------------------------------------#

# Install cellranger-atac-2.1.0.tar.gz
tar -xzvf cellranger-9.0.1.tar.gz

# Install Human reference data
tar -xzvf refdata-gex-GRCh38-2024-A.tar.gz


# add current path to PATH
export PATH=$PATH:/projects/compbio/users/xran2/tool/scRNA_map/cellranger-9.0.1

cellranger sitecheck > sitecheck.txt

cellranger testrun --id=check_install
