// Trying different implementations of popcount
// Can use hardware instruction based on the architecture
// Source: Stack Overflow, Wiki, etc

// the below two functions can be templated
uint32_t popcount32(uint32_t i)
{
  i = i - ((i >> 1) & 0x55555555);
  i = (i & 0x33333333) + ((i >> 2) & 0x33333333);
  return (((i + (i >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;
}

uint32_t popcount64(uint64_t i)
{
  i = i - ((i >> 1) & 0x5555555555555555);
  i = (i & 0x3333333333333333) + ((i >> 2) & 0x3333333333333333);
  i = (i + (i >> 4)) & 0x0f0f0f0f0f0f0f0f;
  i = i + (i >> 8);
  i = i + (i >> 16);
  i = i + (i >> 32);
  return (uint32_t)i & 0x7f;
}

/*
static unsigned char wordbits[65536] = { bitcounts of ints between 0 and 65535 };
static int popcount( unsigned int i )
{
  return( wordbits[i&0xFFFF] + wordbits[i>>16] );
}
*/
