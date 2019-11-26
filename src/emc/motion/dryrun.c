/*
Copyright: 2019
Author: Dewey Garrett <dgarrett@panix.com>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

// debugging using stdio (use with rtpreempt only):
#define DRYRUN_DBG
#undef  DRYRUN_DBG
#ifdef  DRYRUN_DBG //{
#include <stdio.h>
#endif             //}

#include "hal.h"
#include "mot_priv.h"
#include "motion_debug.h"
#include "motion_types.h"
#include <rtapi_math.h>
#include <rtapi_string.h>
#include "dryrun.h"

// strings for translation
#define _(s) (s)

static EmcPose save_pose;
static double  save_joints[EMCMOT_MAX_JOINTS];

static bool       dryrun_saved = 0;
static bool last_dryrun_active = 0;
static bool     dryrun_started = 0;
static bool    dryrun_stopping = 0;
static bool        dryrun_fini = 0;
static int           dryrun_ct = 0;

// hal_data for dryrun outputs:
dryrun_motstatus_data_t  *dryrun_motstatus_data = 0;

//forward declarations
static int  dryrun_export_motstatus(int id, int n);
static void dryrun_once(void);
static void dryrun_clear_stat(void);
static void dryrun_clear_dryrun_stat();

// local vars & functions -------------------------------------
static double *pcmd_p[EMCMOT_MAX_AXIS];
static void dryrun_once(void)
{
    // *pcmd_p[0] is emcmotStatus->carte_pos_cmd.tran.x
    // *pcmd_p[1] is emcmotStatus->carte_pos_cmd.tran.y
    //  etc.
    pcmd_p[0] = &(emcmotStatus->carte_pos_cmd.tran.x);
    pcmd_p[1] = &(emcmotStatus->carte_pos_cmd.tran.y);
    pcmd_p[2] = &(emcmotStatus->carte_pos_cmd.tran.z);
    pcmd_p[3] = &(emcmotStatus->carte_pos_cmd.a);
    pcmd_p[4] = &(emcmotStatus->carte_pos_cmd.b);
    pcmd_p[5] = &(emcmotStatus->carte_pos_cmd.c);
    pcmd_p[6] = &(emcmotStatus->carte_pos_cmd.u);
    pcmd_p[7] = &(emcmotStatus->carte_pos_cmd.v);
    pcmd_p[8] = &(emcmotStatus->carte_pos_cmd.w);
} // dryrun_once()

int dryrun_setup() {
    static int dryrun_comp_id;
    int n,r;

    // hal data for dryrun substitue output pins
    dryrun_motstatus_data = hal_malloc(sizeof(dryrun_motstatus_data_t));
    dryrun_comp_id = hal_init("dryrun");
    if (dryrun_comp_id < 0) {
        rtapi_print_msg(RTAPI_MSG_ERR, _("DRYRUN: hal_init() failed\n"));
        return -1;
    }

    r = dryrun_export_motstatus(dryrun_comp_id, n);
    if (r != 0){
        rtapi_print_msg(RTAPI_MSG_ERR, _("DRYRUN: motstatus export failed\n"));
        return -1;
    }
    hal_ready(dryrun_comp_id);

    return 0; //ok
} // dryrun_setup()

void dryrun_save(void)
{
    KINEMATICS_FORWARD_FLAGS fflags = 0;
    KINEMATICS_INVERSE_FLAGS iflags = 0;
    int joint_num;

    for (joint_num = 0; joint_num < 3+ALL_JOINTS; joint_num++) {
        save_joints[joint_num] = (&joints[joint_num])->pos_fb;
    }

    kinematicsForward(save_joints, &save_pose, &fflags, &iflags);
    dryrun_saved = 1;

    dryrun_show(1,"SAVE");
} // dryrun_save()

bool dryrun_end()
{
    return(dryrun_fini);
} // dryrun_end()

static void dryrun_clear_stat()
{
    *(emcmot_hal_data->coord_mode)     = 0;
    *(emcmot_hal_data->teleop_mode)    = 0;
    *(emcmot_hal_data->in_position)    = 0;
    *(emcmot_hal_data->tp_reverse)     = 0;
    return;
} // dryrun_clear_stat()

static void dryrun_clear_dryrun_stat()
{
    *(dryrun_motstatus_data->coord_mode)  = 0;
    *(dryrun_motstatus_data->teleop_mode) = 0;
    *(dryrun_motstatus_data->in_position) = 0;
    *(dryrun_motstatus_data->tp_reverse)  = 0;
    return;
} // dryrun_clear_dryrun_stat()

void dryrun_manage(void)
{
    static bool gave_disallow_msg = 0;
    static int               once = 1;
    char*  msg = "";
    char*  msg_extra = "";
    if (once) {dryrun_once(); once=0;}
    dryrun_ct++; // for debugging

    // check status from prior cycle
    // if dryrun_fini, done
    if (dryrun_fini) {
        dryrun_show(1,"FINI");
        dryrun_fini = 0;
        if (emcmotExternalOffsetsEnabled() ) {
            reportError("Force machine off"
                        "\nDryrun ended with external_offsets enabled");
            SET_MOTION_ENABLE_FLAG(0);
            emcmotDebug->enabling = 0;
        }
        return;
    }

    // check status from prior cycle
    // if stopping, initiate dryrun_fini cycle
    if (dryrun_stopping) {
        dryrun_show(1,"STOPPING");
        dryrun_stopping = 0;
        dryrun_fini     = 1;
        dryrun_clear_dryrun_stat();
        return;
    };

    // Require motion enable
    if (   !GET_MOTION_ENABLE_FLAG() ) {
        *(emcmot_hal_data->is_dryrun)    = 0;
        *(emcmot_hal_data->isnot_dryrun) = !*(emcmot_hal_data->is_dryrun);
        dryrun_started = 0;
        return;
    }

    // Specialcase, unhomed after starting
    if (  !checkAllHomed()
        && dryrun_started) {
        *(emcmot_hal_data->is_dryrun)    = 0;
        *(emcmot_hal_data->isnot_dryrun) = !*(emcmot_hal_data->is_dryrun);
        dryrun_started = 0;
        rtapi_print_msg(RTAPI_MSG_ERR,"Machine unhomed, dryrun stopped\n");
        return;
    }

    last_dryrun_active = dryrun_active();

    // Disallow startup if motion.dryrun-inhibit-code
    if (   *emcmot_hal_data->dryrun_start
        && *emcmot_hal_data->dryrun_inhibit_code ) {
        msg ="dryrun start prevented by motion.dryrun-inhibit-code\n";
        msg_extra = "";
        switch (*emcmot_hal_data->dryrun_inhibit_code) {
           case DRYRUN_INHIBIT_TOOL_LOADED:
                msg_extra="Inhibit: Tool Loaded"
                          "\nRemove tool first (T0M6)\n";break;
           case DRYRUN_INHIBIT_COOLANT_ON:
                msg_extra="Inhibit: Coolant active\n";   break;
           case DRYRUN_INHIBIT_OTHER:
                msg_extra="Inhibit: Other\n";            break;
           default:
                msg_extra="Inhibit: Unknown code\n";
        } 
        goto disallow_start;
    }

    // Disallow dryrun startup while in coord (running prog/mdi)
    if (   *emcmot_hal_data->dryrun_start
        && !dryrun_started
        && GET_MOTION_COORD_FLAG()) {
        msg = "Cannot start dryrun while mdi or prog"
              "\nHint: 1) If program running, stop it first"
              "\nHint: 2) Force manual mode";
        goto disallow_start;
    }

    // Disallow startup if not homed
    if (   *emcmot_hal_data->dryrun_start
        && !checkAllHomed() ) {
        msg = "Cannot start dryrun until homed";
        goto disallow_start;
    }

    // Disallow startup if external_offsets enabled for any axis
    if (   *emcmot_hal_data->dryrun_start
        && emcmotExternalOffsetsEnabled() ) {
        msg ="Cannot start dryrun with external offsets enabled";
        goto disallow_start;
    }

    // Disallow startup if any spindle is on
    if (   *emcmot_hal_data->dryrun_start
        && emcmotSpindleIsOn() ) {
        msg ="Cannot start dryrun with spindle on ";
        goto disallow_start;
    }

    // Detect disallowed ending of request
    if (   *emcmot_hal_data->dryrun_stop
        && dryrun_started
        && GET_MOTION_COORD_FLAG() ) {
        // cannot distinguish cases in motmod, so give hints:
        msg = "Cannot stop dryrun while mdi or prog"
              "\nHint: 1) If running program, stop it first"
              "\nHint: 2) Force manual mode";
        goto disallow_stop;
    }

    // Disallow stop if spindle on
    if (   *emcmot_hal_data->dryrun_stop
        && dryrun_started
        && emcmotSpindleIsOn() ) {
        // cannot distinguish cases in motmod, so give hints:
        msg = "Cannot stop dryrun with spindle on";
        goto disallow_stop;
    }

    // Allowed startup of dryrun
    if (   *emcmot_hal_data->dryrun_start
        && !dryrun_started) {
       dryrun_clear_stat();
       dryrun_save();
    }

    gave_disallow_msg = 0;
    // Detect allowed ending of request
    if (   *emcmot_hal_data->dryrun_stop
        && dryrun_started) {
       dryrun_stopping = 1;
       SET_MOTION_COORD_FLAG(0);
       SET_MOTION_TELEOP_FLAG(1);
       dryrun_show(1,"STOP_REQUEST");
    }

    if ( *emcmot_hal_data->dryrun_start) { dryrun_started = 1; }
    if ( *emcmot_hal_data->dryrun_stop)  { dryrun_started = 0; }

    *(emcmot_hal_data->is_dryrun)    = dryrun_started;
    *(emcmot_hal_data->isnot_dryrun) = !*(emcmot_hal_data->is_dryrun);
    return;

disallow_start:
    *(emcmot_hal_data->is_dryrun)    = 0;
    *(emcmot_hal_data->isnot_dryrun) = !*(emcmot_hal_data->is_dryrun);
    dryrun_started = 0;
    if (!gave_disallow_msg) {
        rtapi_print_msg(RTAPI_MSG_ERR,msg);
        if (strcmp(msg_extra,"")) {
            rtapi_print_msg(RTAPI_MSG_ERR,msg_extra);
        }
    }
    gave_disallow_msg = 1;
    return;

disallow_stop:
    if (!gave_disallow_msg) {
        rtapi_print_msg(RTAPI_MSG_ERR, msg);
    }
    gave_disallow_msg = 1;
    return;

} // dryrun_manage()

bool dryrun_active(void)
{
    return(   dryrun_started
           || dryrun_stopping
           || dryrun_fini
          );
} // dryrun_active()

void dryrun_restore(void)
{
    double restore_joints[EMCMOT_MAX_JOINTS];
    int joint_num,axis_num;

    if (dryrun_saved) {
        emcmotStatus->carte_pos_cmd = save_pose;
        kinematicsInverse(&emcmotStatus->carte_pos_cmd,restore_joints,&iflags,&fflags);
        for (joint_num = 0; joint_num < ALL_JOINTS; joint_num++) {
            (&joints[joint_num])->pos_cmd = restore_joints[joint_num];
            (&joints[joint_num])->pos_fb  = restore_joints[joint_num];
            // drain coord mode interpolators:
            cubicDrain( &(&joints[joint_num])->cubic);
        }
        for (axis_num = 0; axis_num < EMCMOT_MAX_AXIS; axis_num++) {
            (&axes[axis_num])->teleop_tp.curr_pos = *pcmd_p[axis_num];
            (&axes[axis_num])->teleop_tp.pos_cmd  = *pcmd_p[axis_num];
        }

        dryrun_show(1,"RESTORE");
    } else {
        dryrun_show(1,"RESTORE:UNEXPECTED");
    }
} //dryrun_restore()

void dryrun_show(bool force,char*msg)
{
#ifdef DRYRUN_DBG //{
    char * mstate;
    mstate="UNKNOWN";

    // The 'force' parameter is available for local tests here:
    // Example:
    // if (   !force
    //     && !(testhere) {
    //    return;
    // }

    switch (emcmotStatus->motion_state) {
      case EMCMOT_MOTION_FREE:     mstate="FREE";    break;
      case EMCMOT_MOTION_TELEOP:   mstate="TELEOP";  break;
      case EMCMOT_MOTION_COORD:    mstate="COORD";   break;
      case EMCMOT_MOTION_DISABLED: mstate="DISABLED";break;
    }

    fprintf(stderr,"\n%1s %8d %-20s active=%d stopping=%d fini=%d mstate=%s\n"
           ,force?"F":"*"
           ,dryrun_ct,msg,dryrun_active(),dryrun_stopping,dryrun_fini
           ,mstate
           );
    fprintf(stderr,"simp: %8.4f %8.4f %8.4f cur:   %8.4f %8.4f %8.4f\n"
           ,(&axes[0])->teleop_tp.pos_cmd
           ,(&axes[1])->teleop_tp.pos_cmd
           ,(&axes[2])->teleop_tp.pos_cmd
           ,(&axes[0])->teleop_tp.curr_pos
           ,(&axes[1])->teleop_tp.curr_pos
           ,(&axes[2])->teleop_tp.curr_pos
           );
    fprintf(stderr,"abc:  %8.4f %8.4f %8.4f\n"
           ,emcmotStatus->carte_pos_cmd.a
           ,emcmotStatus->carte_pos_cmd.b
           ,emcmotStatus->carte_pos_cmd.c
           );
    fprintf(stderr,"Pose: %8.4f %8.4f %8.4f\n"
           ,emcmotStatus->carte_pos_cmd.tran.x
           ,emcmotStatus->carte_pos_cmd.tran.y
           ,emcmotStatus->carte_pos_cmd.tran.z
           );
    fprintf(stderr,"cmd:  %8.4f %8.4f %8.4f motor: %8.4f %8.4f %8.4f\n"
           ,(&joints[0])->pos_cmd
           ,(&joints[1])->pos_cmd
           ,(&joints[2])->pos_cmd
           ,(&joints[0])->motor_pos_cmd
           ,(&joints[1])->motor_pos_cmd
           ,(&joints[2])->motor_pos_cmd
           );
    fprintf(stderr,"fb:   %8.4f %8.4f %8.4f fb:    %8.4f %8.4f %8.4f\n"
           ,(&joints[0])->pos_fb
           ,(&joints[1])->pos_fb
           ,(&joints[2])->pos_fb
           ,(&joints[0])->motor_pos_fb
           ,(&joints[1])->motor_pos_fb
           ,(&joints[2])->motor_pos_fb
           );
#endif //} DRYRUN_DBG
    return;
} // dryrun_show()



static int dryrun_export_motstatus(int id, int n)
{
    int r;
    dryrun_motstatus_data_t* p = dryrun_motstatus_data;

    if ((r = hal_pin_bit_newf(  HAL_OUT, &(p->coord_mode),     id, "dmotion.coord-mode")    ) != 0) return r;
    if ((r = hal_pin_bit_newf(  HAL_OUT, &(p->teleop_mode),    id, "dmotion.teleop-mode")   ) != 0) return r;
    if ((r = hal_pin_s32_newf(  HAL_OUT, &(p->motion_type),    id, "dmotion.motion-type")   ) != 0) return r;
    if ((r = hal_pin_bit_newf(  HAL_OUT, &(p->tp_reverse),     id, "dmotion.tp-reverse")    ) != 0) return r;
    if ((r = hal_pin_s32_newf(  HAL_OUT, &(p->program_line),   id, "dmotion.program-line")  ) != 0) return r;
    if ((r = hal_pin_bit_newf(  HAL_OUT, &(p->in_position),    id, "dmotion.in-position")   ) != 0) return r;
    if ((r = hal_pin_float_newf(HAL_OUT, &(p->distance_to_go), id, "dmotion.distance-to-go")) != 0) return r;


    return 0;
} // dryrun_export_motstatus()
