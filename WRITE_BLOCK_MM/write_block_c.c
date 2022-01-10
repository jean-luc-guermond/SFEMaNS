#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <zlib.h>

typedef enum
{
  step_A,
  step_B,
  step_C
} base64_encodestep;

typedef struct
{
  base64_encodestep step;
  char              result;
} base64_encodestate;

void
base64_init_encodestate(base64_encodestate *state_in)
{
  state_in->step   = step_A;
  state_in->result = 0;
}


static inline char
base64_encode_value(char value_in)
{
  static const char *encoding =
    "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

  if (value_in > 63)
    return '=';
  return encoding[(int)value_in];
}

int
base64_encode_block(const char *        plaintext_in,
                    int                 length_in,
                    char *              code_out,
                    base64_encodestate *state_in)
{
  const char *      plainchar    = plaintext_in;
  const char *const plaintextend = plaintext_in + length_in;
  char *            codechar     = code_out;
  char              result;

  result = state_in->result;

  switch (state_in->step)
    {
      while (true)
        {
          case step_A:
            {
              if (plainchar == plaintextend)
                {
                  state_in->result = result;
                  state_in->step   = step_A;
                  return codechar - code_out;
                }
              const char fragment = *plainchar++;
              result              = (fragment & 0x0fc) >> 2;
              *codechar++         = base64_encode_value(result);
              result              = (fragment & 0x003) << 4;
            }
          case step_B:
            {
              if (plainchar == plaintextend)
                {
                  state_in->result = result;
                  state_in->step   = step_B;
                  return codechar - code_out;
                }
              const char fragment = *plainchar++;
              result |= (fragment & 0x0f0) >> 4;
              *codechar++ = base64_encode_value(result);
              result      = (fragment & 0x00f) << 2;
            }
          case step_C:
            {
              if (plainchar == plaintextend)
                {
                  state_in->result = result;
                  state_in->step   = step_C;
                  return codechar - code_out;
                }
              const char fragment = *plainchar++;
              result |= (fragment & 0x0c0) >> 6;
              *codechar++ = base64_encode_value(result);
              result      = (fragment & 0x03f) >> 0;
              *codechar++ = base64_encode_value(result);
            }
        }
    }
  /* control should not reach here */
  return codechar - code_out;
}

int
base64_encode_blockend(char *code_out, base64_encodestate *state_in)
{
  char *codechar = code_out;

  switch (state_in->step)
    {
      case step_B:
        *codechar++ = base64_encode_value(state_in->result);
        *codechar++ = '=';
        *codechar++ = '=';
        break;
      case step_C:
        *codechar++ = base64_encode_value(state_in->result);
        *codechar++ = '=';
        break;
      case step_A:
        break;
    }
  *codechar++ = '\0';

  return codechar - code_out;
}

/*
 * We use the bb64 library for base64 encoding:
 */

char *
encode_block(const char *data, const int data_size)
{
  base64_encodestate state;
  base64_init_encodestate(&state);

  char *encoded_data = malloc(2 * data_size + 1);

  const int encoded_length_data =
    base64_encode_block(data, data_size, encoded_data, &state);
  base64_encode_blockend(encoded_data + encoded_length_data, &state);

  return encoded_data;
}

/*
 * We use the zlib library for compression:
 */

void
write_compressed_block(int *fd, void **block, int *block_size)
{
  // Compress data:

  uLongf compressed_data_length = compressBound(*block_size);
  char * compressed_data        = malloc(compressed_data_length);

  compress2((Bytef *)compressed_data,
            &compressed_data_length,
            (const Bytef *)*block,
            *block_size,
            Z_BEST_SPEED);
/*             Z_NO_COMPRESSION); */

  const uint32_t compression_header[4] = {
    1,                                 /* number of blocks */
    (uint32_t)*block_size,             /* size of block */
    (uint32_t)*block_size,             /* size of last block */
    (uint32_t)compressed_data_length}; /* list of compressed sizes of blocks */

  char *encoded_header =
    encode_block((const char *)(&compression_header[0]),
                 4 * sizeof(compression_header[0]));

  write(*fd, encoded_header, strlen(encoded_header));
  free(encoded_header);

  char *encoded_data = encode_block(compressed_data, compressed_data_length);
  free(compressed_data);

  write(*fd, encoded_data, strlen(encoded_data));
  free(encoded_data);

  write(*fd, "\n", 1);
}
