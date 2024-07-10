        character*(1) function MYCHAR1(i)
        integer i
         MYCHAR1(1:1) = ACHAR(ICHAR('0')+i)
        end function

        character*(2) function MYCHAR2(i)
        integer i
         MYCHAR2(1:1) = ACHAR(ICHAR('0')+INT(i/10))
         MYCHAR2(2:2) = ACHAR(ICHAR('0')+i-10*INT(i/10))
        end function

        character*(3) function MYCHAR3(i)
        integer i
         MYCHAR3(1:1) = ACHAR(ICHAR('0')+INT(i/100))
         MYCHAR3(2:2) = ACHAR(ICHAR('0')+INT(i/10)-10*INT(i/100))
         MYCHAR3(3:3) = ACHAR(ICHAR('0')+i-10*INT(i/10))
        end function

        character*(4) function MYCHAR4(i)
        integer i
         MYCHAR4(1:1) = ACHAR(ICHAR('0')+INT(i/1000))
         MYCHAR4(2:2) = ACHAR(ICHAR('0')+INT(i/100)-10*INT(i/1000))
         MYCHAR4(3:3) = ACHAR(ICHAR('0')+INT(i/10)-10*INT(i/100))
         MYCHAR4(4:4) = ACHAR(ICHAR('0')+i-10*INT(i/10))
        end function

        character*(5) function MYCHAR5(i)
        integer i
         MYCHAR5(1:1) = ACHAR(ICHAR('0')+INT(i/10000))
         MYCHAR5(2:2) = ACHAR(ICHAR('0')+INT(i/1000)-10*INT(i/10000))
         MYCHAR5(3:3) = ACHAR(ICHAR('0')+INT(i/100)-10*INT(i/1000))
         MYCHAR5(4:4) = ACHAR(ICHAR('0')+INT(i/10)-10*INT(i/100))
         MYCHAR5(5:5) = ACHAR(ICHAR('0')+i-10*INT(i/10))
        end function

        character*(6) function MYCHAR6(i)
        integer i
         MYCHAR6(1:1) = ACHAR(ICHAR('0')+INT(i/100000))
         MYCHAR6(2:2) = ACHAR(ICHAR('0')+INT(i/10000)-10*INT(i/100000))
         MYCHAR6(3:3) = ACHAR(ICHAR('0')+INT(i/1000)-10*INT(i/10000))
         MYCHAR6(4:4) = ACHAR(ICHAR('0')+INT(i/100)-10*INT(i/1000))
         MYCHAR6(5:5) = ACHAR(ICHAR('0')+INT(i/10)-10*INT(i/100))
         MYCHAR6(6:6) = ACHAR(ICHAR('0')+i-10*INT(i/10))
        end function

