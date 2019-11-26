#ifndef DRYRUN_H
#define DRYRUN_H
#include <rtapi_bool.h>

//-------------------------------------------------
// Dryrun messages (dryrun_helper handles)
// mask out using pin: dryrun_helper.message-mask

// dryrun_helper internal messages:
#define DRYRUN_START_MSG           0x01
#define DRYRUN_STOP_MSG            0x02

// Dryrun io messages from iocontrol to dryrun_helper
//   From: <= iocontrol.0.dryrun-io-message
//   To:   => dryrun_helper.io-message
#define DRYRUN_TOOL_PREPARED_MSG   0x04
#define DRYRUN_TOOL_CHANGED_MSG    0x08
#define DRYRUN_TOOL_REMAINING_MSG  0x10
#define DRYRUN_LAST_MSG_NUMBER     0x10

//-------------------------------------------------
// Inhibit codes: disallow dryrun startup
//  From: <= dryrun_helper.dryrun-inhibit-code
//  To:   => motion.dryrun-inhibit-code
#define DRYRUN_INHIBIT_COOLANT_ON  0x01
#define DRYRUN_INHIBIT_TOOL_LOADED 0x02
#define DRYRUN_INHIBIT_OTHER       0x80
//-------------------------------------------------

// dryrun speed warping:
#define DRYRUN_MAX_RUNCYCLES 100

int dryrun_setup(void);    // onetime

void dryrun_manage(void);  // per cycle
bool dryrun_end(void);     // query
bool dryrun_active(void);  // query
void dryrun_restore(void); // restore state

// for debugging (for rtpreempt, uses stdio):
void dryrun_show(bool force,char* msg);

//-----------------------------------------------------
// support for (some) dryrun motstatus out pins (dmotion.*)
typedef struct {
    hal_bit_t   *coord_mode;
    hal_bit_t   *teleop_mode;
    hal_s32_t   *motion_type;
    hal_bit_t   *tp_reverse;
    hal_s32_t   *program_line;
    hal_bit_t   *in_position;
    hal_float_t *distance_to_go;
} dryrun_motstatus_data_t;

extern dryrun_motstatus_data_t* dryrun_motstatus_data;

#endif /* DRYRUN_H */
