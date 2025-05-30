//------------------------------------------------------------------------------
// This file is part of the OpenStructure project <www.openstructure.org>
//
// Copyright (C) 2008-2020 by the OpenStructure authors
//
// This library is free software; you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License as published by the Free
// Software Foundation; either version 3.0 of the License, or (at your option)
// any later version.
// This library is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
// details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library; if not, write to the Free Software Foundation, Inc.,
// 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
//------------------------------------------------------------------------------

/*
  Author: Stefan Scheuber
 */



#include "seq_selection_painter.hh"
#include <QtGui>

namespace ost { namespace gui {

SeqSelectionPainter::SeqSelectionPainter(QObject* parent)
    : Painter(parent), focus_color_(255,0,0,0), mouse_over_color_(240,240,240,192)
{

}

void SeqSelectionPainter::Paint(QPainter* painter, const QStyleOptionViewItem& option, const QModelIndex& index){
  painter->save();
  /*
  if (option.state & QStyle::State_HasFocus){
    painter->fillRect(option.rect, focus_color_);
  }
  */
  if (option.state & QStyle::State_MouseOver){
    painter->fillRect(option.rect, mouse_over_color_);
  }
  if (option.state & QStyle::State_Selected){
    QColor color = option.palette.highlight().color();
    color.setAlpha(128);
    painter->fillRect(option.rect, color);
  }
  painter->restore();
}

}}
