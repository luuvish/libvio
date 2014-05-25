{
  "targets": [
    {
      "target_name": "codec",
      "include_dirs": [
        "h264/core",
        "h264/decoder",
        "h264/framebuf",
        "h264/parser"
      ],
      "sources": [
        "codec.cc",
        "h264/core/input_parameters.cc",
        "h264/core/ldecod.cc",
        "h264/core/report.cc",
        "h264/core/slice_data.cc",
        "h264/core/slice_fmo.cc",
        "h264/core/slice_header.cc",
        "h264/core/slice_ref_list.cc",
        "h264/decoder/deblock.cc",
        "h264/decoder/decoder.cc",
        "h264/decoder/erc_api.cc",
        "h264/decoder/inter_prediction.cc",
        "h264/decoder/intra_prediction.cc",
        "h264/decoder/transform.cc",
        "h264/framebuf/dpb.cc",
        "h264/framebuf/dpb_erc.cc",
        "h264/framebuf/image_data.cc",
        "h264/framebuf/memalloc.cc",
        "h264/framebuf/output.cc",
        "h264/framebuf/picture.cc",
        "h264/parser/bitstream.cc",
        "h264/parser/bitstream_cabac.cc",
        "h264/parser/bitstream_rtp.cc",
        "h264/parser/interpret.cc",
        "h264/parser/interpret_mb.cc",
        "h264/parser/interpret_mv.cc",
        "h264/parser/interpret_rbsp.cc",
        "h264/parser/interpret_residual.cc",
        "h264/parser/interpret_se.cc",
        "h264/parser/interpret_sei.cc",
        "h264/parser/neighbour.cc"
      ],
      "conditions": [
        [ 'OS=="mac"', {
          "xcode_settings": {
            "OTHER_CPLUSPLUSFLAGS" : ["-std=c++11", "-stdlib=libc++"],
            "OTHER_LDFLAGS": ["-stdlib=libc++"],
            "MACOSX_DEPLOYMENT_TARGET": "10.7"
          }
        }]
      ]
    }
  ]
}
