BEGIN { print "Don’t Panic!" }
# procura por foo em data
awk '/foo/ { print $0 }' data

# acha a maior linha
awk '{ if (length($0) > max) max = length($0) } END { print max }' data

# = anterior usando expand
expand data | awk '{ if (x < length()) x = length() } END { print "maximum line length is " x }'

# pega as linhas maiores que 80
awk 'length($0) > 80' data

# pega as linhas com pelo menos um campo
awk 'NF > 0' data

# Gera números aleatórios de 0 a 100 inclusive
awk 'BEGIN { for (i = 1; i <= 7; i++) print int(101 * rand()) }'

# usa o pipe para capturar a 5 coluna e inprimir
ls -l files | awk '{ x += $5 } END { print "total K-bytes: " (x + 1023)/1024 }'

# Pega a primeira coluna separada por : e mostra pelo sort
awk -F: '{ print $1 }' /etc/passwd | sort

# numero de linhas no arquivo
awk 'END { print NR }' data

# pega as linhas impares (para pares == 1)
awk 'NR % 2 == 0' data

# se a coluna 6 for Nov soma na coluna 5
ls -l | awk '$6 == "Nov" { sum += $5 } END { print sum }'

# usa \ para continuar (quando dentro de script
awk '/This regular expression is too long, so continue it\
 on the next line/ { print $1 }'

# pega todas as linhas da 6 coluna que tem J
ls -l | awk '{ if ($6 ~ /J/) print }'

# pega todas as linhas da 6 coluna que não tem J
ls -l | awk '{ if ($6 !~ /J/) print }'

# usa \ para caracteres especiais como " e \ (\a produz um som)
awk 'BEGIN { print "He said \"hi!\" to the bar \\." }'

# pega as colunas 1 e 3 do arquivos
awk '{print $1,$3}' my.dat

# conta palavras ou campos por linha
awk '{print NF}' my.dat

# mostra formatado notação em C
awk '{printf "%4d %8.4f\n",  $1, $2}' my.dat

# suma para a coluna 1 
awk '{s+=$1} END {print s}' my.dat

# comandos awk em um arquivo
awk -f awkfile filename

# mesmo que 666
awk '/666/'  filename

# linhas com alphanumeric
awk '/[[:alnum:]]/'  filename

# linhas com ...999... na coluna 2
awk '$2 ~/999/'  filename

# linhas com exatamente 999 na coluna 2
awk '$2 =="999"'  filename

# mostra cada 10 linhas 
awk 'print (NR%10)' filename

# pega o 10th campo quando alguns campos tem valores perdidos ou sem espaço entre eles
awk 'BEGIN {FIELDWIDTHS="5 5 5 5 2 3 2 2 9 8 8 8 8 8"} {print $10}' filename