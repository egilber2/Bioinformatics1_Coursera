
def PatternCount(Text, Pattern):
    """Returns count of Pattern occurrences in a string Text.
    
    Keyword arguments:
    ------------------
    Text -- string of DNA
    Pattern -- pattern of DNA of interest
    
    Returns:
    --------
    int
    """
    count=0
    j=len(Pattern)
    for i in range(len(Text)-len(Pattern)+1):
        if Text[i: i+j]==Pattern:
            count+=1
        else:
            i+=1
    return count


def FrequentWords(Text, k):
    """Returns the most frequent k-mer in a string.
    
    Keyword arguments:
    ------------------
    Text -- the string of DNA
    k -- an integer for length of k-mer
    
    Returns:
    --------
    list   
    """
    
    FrequentPatterns = []
    count_list=[1]
    maxCount=0
    for i in range(len(Text)-k +1):
        Pattern = Text[i:i + k]
        count = PatternCount(Text, Pattern)
        if count >= max(count_list):
            count_list.append(count)
            FrequentPatterns = []
            FrequentPatterns.append(Pattern)
        elif count == max(count_list):
            FrequentPatterns.append(Pattern)
        else:
            continue
    a = list(set(FrequentPatterns))
    return a 



def revComplement(p):
    """Returns the reverse complement for strand of DNA.
    
    Keyword arguments:
    ------------------
    p -- string of DNA
    
    Returns:
    --------
    string 
    """
    revComp=[]
    revString = p[::-1]
    for char in revString:
        if char=='A':
            revComp.append('T')
        elif char=='T':
            revComp.append('A')
        elif char=='C':
            revComp.append('G')
        elif char=='G':
            revComp.append('C')
    a = "".join(revComp)
    return a


def patternMatch(pattern, genome):
    """Returns index where pattern starts in genome.
    
    Keyword arguments:
    ------------------
    pattern -- string, DNA pattern
    genome -- string, DNA
    
    Returns:
    --------
    Position index of genome where pattern begins
    """
    indexList=[]
    j=len(pattern)
    for i in range(len(genome)-len(pattern)+1):
        if genome[i: i+j]==pattern:
            indexList.append(i)
        else:
            continue
    a = print(*indexList, sep=' ')
    return a


def FrequentWordsII(Text, k):
    FrequentPatterns = []
    count_list=[1]
    maxCount=0
    for i in range(len(Text)-k +1):
        Pattern = Text[i:i + k]
        count = PatternCount(Text, Pattern)
        if count >= max(count_list):
            count_list.append(count)
            FrequentPatterns = []
            FrequentPatterns.append(Pattern)
        elif count == max(count_list):
            FrequentPatterns.append(Pattern)
        else:
            continue
    #x = list(set(FrequentPatterns))
    return  FrequentPatterns


def clump_find(genome, k, L, t):
    clump_list=[]
    for i in range(len(genome)-L +1):
        Text=genome[i:i+L]
        clump = FrequentWordsII(Text, k)
        z = clump_list + clump
    #z = list(set(clump_list))
    return z
        

def symbolToNumber(symbol):
    """Converts DNA base to corresponding number.
    
    Keyword arguments:
    ------------------
    symbol -- DNA base
    
    Returns:
    --------
    int
    """
    if symbol=='A':
        x = 0
    elif symbol=='C':
        x= 1
    elif symbol == 'G':
        x = 2
    elif symbol == 'T':
        x = 3
    return x



def NumberToSymbol(num):
    """Coverts int to corresponding DNA base abbreviation.
    
    Keyword arguments:
    ------------------
    num -- int 0,1,2, or 3
    
    Returns:
    --------
    string -- DNA base abbreviation
    """
    if num == 0:
        x = 'A'
    elif num == 1:
        x = 'C'
    elif num == 2:
        x = 'G'
    elif num == 3:
        x = 'T'
    return x


def PatternToNumber(pattern):
    '''Converts DNA string to base 4 number.
    
    Keyword arguments:
    ------------------
    pattern -- string, DNA
    
    Returns:
    --------
    int -- base 4 number corresponding to DNA string    
    '''
    patternList = list(pattern)
    if not patternList:
        return 0
    symbol=pattern[-1]
    prefix = pattern[:-1]
    return 4*PatternToNumber(prefix) + symbolToNumber(symbol)
    

def NumberToPattern(index, k):
    """Converts an int to corresponding string of DNA.
    
    Keyword arguments:
    ------------------
    index -- int
    k -- int, length of k-mer
    
    Returns:
    --------
    string of DNA bases

    """
    if k == 1:
        return NumberToSymbol(index)
    prefixIndex = index//4
    r = index % 4
    symbol = NumberToSymbol(r)
    PrefixPattern = NumberToPattern(prefixIndex, k-1)
    return PrefixPattern + symbol
    
    

def ComputingFrequencies(Text, k):
    """Creates frequency array k-mers in DNA string.
    
    Keyword arguments:
    ------------------
    Text -- string of DNA
    k -- int for k-mer length
    
    Returns:
    --------
    int -- frequency array 
    """
    FrequencyArray=[]
    for i in range(4**k):
        FrequencyArray.insert(i,0)
    for i in range(len(Text)-k+1):
        pattern = Text[i:i+k]
        j = PatternToNumber(pattern)
        FrequencyArray[j] += 1
    return FrequencyArray


def FasterFrequentWords(Text, k):
    """Returns the most frequent k-mer in a string.
    
    Keyword arguments:
    ------------------
    Text -- the string of DNA
    k -- an integer for length of k-mer
    
    Returns:
    --------
    list of most frequent k-mers  
    """
    FrequencyPatterns = set()
    FrequencyArray = ComputingFrequencies(Text, k)
    maxCount = max(FrequencyArray)
    for i in range(4**k - 1):
        if FrequencyArray[i] == maxCount:
            pattern = NumberToPattern(i, k)
            FrequencyPatterns.add(pattern)
    return FrequencyPatterns


def FindingFrequentWordsBySorting(Text , k):
    """Returns the most frequent k-mer in a DNA string using a sorted 
    frequency array.
    
    Keyword arguments:
    ------------------
    Text -- the string of DNA
    k -- an integer for length of k-mer
    
    Returns:
    --------
    list of most frequent k-mers
    """
    
    FrequentPatterns = set()
    index = []
    count = []
    for i in range(len(Text)-k):
        index.insert(i,0)
        count.insert(i,0)
    for i in range(len(Text)-k):
        pattern = Text[i : i+k]
        index[i] = PatternToNumber(pattern)
        count[i] = 1
    sortedIndex = sorted(index)
    for i in range(len(Text)-k):
        if sortedIndex[i] == sortedIndex[i-1]:
            count[i] = count[i-1] + 1
    maxCount = max(count)
    for i in range(len(Text)-k):
        if count[i] == maxCount:
            pattern = NumberToPattern(sortedIndex[i], k)
            FrequentPatterns.add(pattern)
    return FrequentPatterns


def ClumpFinding(genome, k, t, L):
    """Identifies (L,t)-clumps that occur at least t times in the genome
    in a window of length L.
    
     Keyword arguments:
    ------------------
    genome -- string of DNA bases
    k -- int for k-mer length
    t -- int for minimum number of times k-mer appears in genome window L
    L -- int for length of window searched in genome
    
    Returns:
    --------
    string -- (L,t)-clumps of k-mers in genome
    """
    FrequentPatterns = set()
    clump = []
    for i in range(4**k - 1):
        clump.insert(i, 0)
    for i in range(len(genome)-L):
        Text = genome[i:i+L]
        FrequencyArray = ComputingFrequencies(Text, k)
        for j in range(4**k - 1):
            if FrequencyArray[j] >= t:
                clump[j] = 1
    for i in range(4**k - 1):
        if clump[i] == 1:
            pattern = NumberToPattern(i, k)
            FrequentPatterns.add(pattern)
    return print(*FrequentPatterns, sep = ' ')
            

        
def BetterClumpFinding(genome, k, t, L):
    """Identifies (L,t)-clumps that occur at least t times in the genome
    in a window of length L utilizing a frequency array.
    
     Keyword arguments:
    ------------------
    genome -- string of DNA bases
    k -- int for k-mer length
    t -- int for minimum number of times k-mer appears in genome window L
    L -- int for length of window searched in genome
    
    Returns:
    --------
    string -- (L,t)-clumps of k-mers in genome
    """
    
    FrequentPatterns = set()
    clump = []
    for i in range(4**k - 1):
        clump.insert(i, 0)
    Text = genome[0:L]
    FrequencyArray = ComputingFrequencies(Text, k)
    for i in range(4**k - 1):
        if FrequencyArray[i] >= t:
            clump[i] = 1
    for i in range(len(genome) - L):
        firstPattern = genome[i-1: i-1+k]
        j = PatternToNumber(firstPattern)
        FrequencyArray[j] = FrequencyArray[j] - 1
        lastPattern = genome[i -1 + L - k: i -1 + L]
        j = PatternToNumber(lastPattern)
        FrequencyArray[j] = FrequencyArray[j] + 1
        if FrequencyArray[j] >= t:
            clump[j] = 1
    for i in range(4**k - 1):
        if clump[i] == 1:
            pattern = NumberToPattern(i, k)
            FrequentPatterns.add(pattern)
    return print(*FrequentPatterns, sep = ' ')

