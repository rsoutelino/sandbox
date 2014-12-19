'sed' is a powerful tool for making the same change to many files.

# change (first if any) occurrence on any line of TD to AB 
sed 's/TD/AB/' my.dat

# change (first if any) occurrence on any line of 1 to 2
sed 's/1/2/' my.dat

# change globally occurrences of 1 to 2
sed 's/1/2/g' my.dat

# change first occurrence on any line of any character from A-L to Z
sed 's/[A-L]/Z/' my.dat

# change any 2 in last position of line to 0
sed 's/2$/0/' my.dat

# For first 2 adjacent upper cases letters, reverse and put a space in between.
sed 's/\([A-Z]\)\([A-Z]\)/\2 \1/' my.dat

# Put space after first digit
sed 's/\([0-9]\)/\1 /' my.dat

# Put space after first letter or digit
sed 's/\([0-9,A-Z,a-z]\)/\1 /' my.dat

# Put space after every digit
sed 's/\([0-9]\)/\1 /g' my.dat

# Duplicate every digit
sed 's/\([0-9]\)/\1\1/g' my.dat

# Duplicate every word
sed 's/\([0-9,a-z,A-Z]*[ ]*\)/\1 \1/g' my.dat

# convert dates from mm/dd/yy notation to yy/mm/dd form (for sorting by dates)
sed 's/\([0-9][0-9]\)\/\([0-9][0-9]\)\/\([0-9][0-9]\)/\3\1\2/g' my.dat

# Ele serve para eliminar os códigos de acentuação que programas visuais como o DreanWeaver usam no código HTML no lugar dos caracteres acentuados propriamente ditos. 
    #!/bin/sed -f
    s/á/á/g
    s/ç/ç/g
    s/ê/ê/g
    s/é/é/g
    s/í/í/g
    s/ã/ã/g
    s/õ/õ/g
    s/ó/ó/g
    s/ú/ú/g
    s/à/à/g
    s/"/"/g