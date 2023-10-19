/**
 * @file earlyloc_report.h
 * @author Benjamin Yang in Department of Geology, National Taiwan University
 * @brief
 * @version 0.1
 * @date 2023-10-05
 *
 * @copyright Copyright (c) 2023
 *
 */
#pragma once
/* Earthworm platform include */
#include <transport.h>
/* */
#include <earlyloc.h>

/* */
void el_report_file_output( const HYPO_STATE *, const char *, const char *, const int );
void el_report_ring_output( const HYPO_STATE *, SHM_INFO *, MSG_LOGO *, const char *, const int );
