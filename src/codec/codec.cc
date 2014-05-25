#include <node.h>

using namespace v8;

#include <stdio.h>
#include <string.h>

#include "input_parameters.h"
#include "h264decoder.h"


static void Configure(InputParameters *p_Inp, int ac, char *av[])
{
    memset(p_Inp, 0, sizeof(InputParameters));
    strcpy(p_Inp->infile, "test.264");
    strcpy(p_Inp->outfile, "test_dec.yuv");
    strcpy(p_Inp->reffile, "test_rec.yuv");

    p_Inp->ParseCommand(ac, av);

    fprintf(stdout, "----------------------------- JM %s %s -----------------------------\n", VERSION, EXT_VERSION);
    if (!p_Inp->bDisplayDecParams) {
        fprintf(stdout, "--------------------------------------------------------------------------\n");
        fprintf(stdout, " Input H.264 bitstream                  : %s \n", p_Inp->infile);
        fprintf(stdout, " Output decoded YUV                     : %s \n", p_Inp->outfile);
        fprintf(stdout, " Input reference file                   : %s \n", p_Inp->reffile);
        fprintf(stdout, "--------------------------------------------------------------------------\n");
    }
}


int start(int argc, char** argv)
{
    int iRet;
    int iFramesDecoded = 0;
    InputParameters InputParams;
    DecoderParams Decoder;

    //get input parameters;
    Configure(&InputParams, argc, argv);
    //open decoder;
    Decoder.OpenDecoder(&InputParams);

    //decoding;
    do {
        iRet = Decoder.DecodeOneFrame();
        if (iRet == DEC_EOS || iRet == DEC_SUCCEED)
            iFramesDecoded++;
        else
            //error handling;
            fprintf(stderr, "Error in decoding process: 0x%x\n", iRet);
    } while ((iRet == DEC_SUCCEED) &&
             (Decoder.p_Inp->iDecFrmNum == 0 || iFramesDecoded < Decoder.p_Inp->iDecFrmNum));

    Decoder.FinitDecoder();
    Decoder.CloseDecoder();

    printf("%d frames are decoded.\n", iFramesDecoded);
    return 0;
}


Handle<Value> Main(const Arguments& args) {
  HandleScope scope;

  if (args.Length() < 2) {
    ThrowException(Exception::TypeError(String::New("Wrong number of arguments")));
    return scope.Close(Undefined());
  }

  if (!args[0]->IsString() || !args[1]->IsString()) {
    ThrowException(Exception::TypeError(String::New("Wrong arguments")));
    return scope.Close(Undefined());
  }

  String::Utf8Value arg0(args[0]);
  String::Utf8Value arg1(args[1]);

  const int argc = 5;
  char* argv[] = {"h264", "-i", *arg0, "-o", *arg1};
  start(argc, argv);

  return scope.Close(Undefined());
}

void Init(Handle<Object> exports, Handle<Object> module) {
  Local<Object> h264 = Object::New();
  h264->Set(String::NewSymbol("main"), FunctionTemplate::New(Main)->GetFunction());
  Local<Object> hevc = Object::New();
  hevc->Set(String::NewSymbol("main"), FunctionTemplate::New(Main)->GetFunction());
  Local<Object> vc1 = Object::New();
  vc1->Set(String::NewSymbol("main"), FunctionTemplate::New(Main)->GetFunction());
  Local<Object> vp8 = Object::New();
  vp8->Set(String::NewSymbol("main"), FunctionTemplate::New(Main)->GetFunction());
  Local<Object> vp9 = Object::New();
  vp9->Set(String::NewSymbol("main"), FunctionTemplate::New(Main)->GetFunction());

  exports->Set(String::NewSymbol("h264"), h264);
  exports->Set(String::NewSymbol("hevc"), hevc);
  exports->Set(String::NewSymbol("vc1"), vc1);
  exports->Set(String::NewSymbol("vp8"), vp8);
  exports->Set(String::NewSymbol("vp9"), vp9);
}

NODE_MODULE(codec, Init)
