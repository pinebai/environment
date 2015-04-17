#!python

def anagram(sentence):
    anagrams={}
    words=sentence.split()
    for x in words:
        sx=sorted(x)
        sx=''.join(sx)
        if sx in anagrams:
            if x not in anagrams[sx]:
                anagrams[sx].append(x)
        else:
            anagrams[sx]=[x]
    for key in anagrams:
        if len(anagrams[key])>1:
            print anagrams[key]

sentence="I may opt for a top yam for amy"
anagram(sentence)
