void swapbytes(void *object, size_t size);
void swap2us(unsigned short x);
void swap4ul(unsigned long x);
float swap4f(float f);
void swap2us(unsigned short x)
   {
       unsigned short out = 0;
       int i;
      for(i = 0; i < 2; ++i)
   {
    const unsigned short byte = (x >> 8 * i) & 0xff;
    out |= byte << (8 - 8 * i);
   }
  x=out;
   }

   void swap4ul(unsigned long x)
   {
       unsigned long out = 0;
    int i;
      for(i = 0; i < 4; ++i)
   {
    const unsigned long byte = (x >> 8 * i) & 0xff;
    out |= byte << (24 - 8 * i);
   }
  x=out;
   }
void swapbytes(void *object, size_t size)
{
   unsigned char *start, *end;
   for ( start = object, end = start + size - 1; start < end; ++start, --end )
   {
      unsigned char swap = *start;
      *start = *end;
      *end = swap;
   }
}
   float swap4f(float f)
{
   float retVal;
   char *floatToConvert = ( char* ) & f;
   char *returnFloat = ( char* ) & retVal;
   returnFloat[0] = floatToConvert[3];
   returnFloat[1] = floatToConvert[2];
   returnFloat[2] = floatToConvert[1];
   returnFloat[3] = floatToConvert[0];
   return retVal;
}



