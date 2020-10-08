/******************************************************************************
 OSRA: Optical Structure Recognition Application

 Created by Igor Filippov, 2007-2013 (igor.v.filippov@gmail.com)

 This program is free software; you can redistribute it and/or modify it under
 the terms of the GNU General Public License as published by the Free Software
 Foundation; either version 2 of the License, or (at your option) any later
 version.

 This program is distributed in the hope that it will be useful, but WITHOUT ANY
 WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
 PARTICULAR PURPOSE.  See the GNU General Public License for more details.

 You should have received a copy of the GNU General Public License along with
 this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
 St, Fifth Floor, Boston, MA 02110-1301, USA
 *****************************************************************************/
// File: osra_segment.cpp
//
// Defines page segmentation functions
//

#include <iostream> // std::ostream, std::cout
#include <fstream> 
#include <queue>

#include "osra.h"
#include "osra_common.h"
#include "osra_segment.h"
#include "osra_labels.h"
#include "osra_ocr.h"


void find_connected_components(const Image &image, double threshold, const ColorGray &bgColor,
                               std::vector<std::list<point_t> > &segments,
                               std::vector<std::vector<point_t> > &margins, bool adaptive)
{
  point_t p;
  std::list<point_t> points;
  int speckle_area = 2;
  if (adaptive)
    {
      int speckle_side = std::min(image.columns(), image.rows()) / 200;
      speckle_area = speckle_side * speckle_side;
      if (speckle_area < 2) speckle_area = 2;
    }

  std::vector<std::vector<int> > tmp(image.columns(), std::vector<int> (image.rows(), 0));

  for (unsigned int i = 0; i < image.columns(); i++)
    for (unsigned int j = 0; j < image.rows(); j++)
      if (get_pixel(image, bgColor, i, j, threshold) == 1) // populate with low threshold for future anisotropic smoothing
        tmp[i][j] = 1;


  for (unsigned int i = 0; i < image.columns(); i++)
    for (unsigned int j = 0; j < image.rows(); j++)
      if (tmp[i][j] == 1)
        {
          tmp[i][j] = 2;
          p.x = i;
          p.y = j;
          points.push_back(p);
          std::list<point_t> new_segment;
          std::vector<point_t> new_margin;
          int counter = 0;
          point_t p1;
          while (!points.empty())
            {
              p = points.back();
              points.pop_back();
              new_segment.push_back(p);
              tmp[p.x][p.y] = -1;
              bool on_the_margin = false;

              // "k" should be in range "[0 .. image.columns) intercepted with [p.x - 1 .. p.x + 2)" ==> "p.x + 2" should be positive ==> "p.x >= -1"
              // "l" should be in range "[0 .. image.rows)    intercepted with [p.y - 1 .. p.y + 2)" ==> "p.y + 2" should be positive ==> "p.y >= -1"
              if (p.x >= -1 && p.y >= -1)
                {
                  unsigned int x_lower = p.x > 1 ? p.x - 1 : 0; // "k" cannot be less then zero
                  unsigned int y_lower = p.y > 1 ? p.y - 1 : 0; // "l" cannot be less then zero
                  unsigned int x_upper = p.x + 2;
                  if (x_upper > image.columns())
                    x_upper = image.columns();
                  unsigned int y_upper = p.y + 2;
                  if (y_upper > image.rows())
                    y_upper = image.rows();

                  for (int k = x_lower; k < x_upper; k++)
                    for (int l = y_lower; l < y_upper; l++)
                      {
                        if (tmp[k][l] == 1)
                          {
                            p1.x = k;
                            p1.y = l;
                            points.push_back(p1);
                            tmp[k][l] = 2;
                          }
                        else if ((int) k != p.x && (int) l != p.y && tmp[k][l] == 0)
                          {
                            on_the_margin = true;
                          }
                      }
                }

              if (on_the_margin && (new_margin.size() < PARTS_IN_MARGIN || (counter % PARTS_IN_MARGIN) == 0))
                new_margin.push_back(p);
              if (on_the_margin)
                counter++;
            }
          if (segments.size() > MAX_SEGMENTS)
            return;
          if (new_segment.size() > speckle_area)
            {
              segments.push_back(new_segment);
              margins.push_back(new_margin);
            }
        }
}

void build_explicit_clusters(
    const std::list<std::list<int> > &clusters,
    const std::vector<std::list<point_t> > &segments,
    std::list<std::list<std::list<point_t> > > &explicit_clusters)
{
  explicit_clusters.clear();
  for (std::list<std::list<int> >::const_iterator c = clusters.begin(); c != clusters.end(); c++)
    {
      std::list<std::list<point_t> > set_of_segments;
      for (std::list<int>::const_iterator s = c->begin(); s != c->end(); s++)
        if (!segments[*s].empty())
          set_of_segments.push_back(segments[*s]);
      if (!set_of_segments.empty())
        explicit_clusters.push_back(set_of_segments);
    }
}

void remove_separators(std::vector<std::list<point_t> > &segments,
                       std::vector<std::vector<point_t> > &margins,
                       double max_aspect, unsigned int size)
{
  std::vector<std::list<point_t> >::iterator s;
  std::vector<std::vector<point_t> >::iterator m;
  s = segments.begin();
  m = margins.begin();

  while (s != segments.end() && m != margins.end())
    {
      if (s->size() <= size)
        {
          s++;
          m++;
          continue;
        }

      int stop = INT_MAX, sleft = INT_MAX, sbottom = 0, sright = 0;
      for (std::list<point_t>::iterator p = s->begin(); p != s->end(); p++)
        {
          if (p->x < sleft)
            sleft = p->x;
          if (p->x > sright)
            sright = p->x;
          if (p->y < stop)
            stop = p->y;
          if (p->y > sbottom)
            sbottom = p->y;
        }
      double aspect = 0;

      if (sright != sleft)
        aspect = 1. * (sbottom - stop+1) / (sright - sleft+1); // where did right and left come from?
      if (aspect > max_aspect || aspect < 1. / max_aspect)
        {
          s = segments.erase(s);
          m = margins.erase(m);
        }
      else
        {
          s++;
          m++;
        }
    }
}

double cos_angle_between_points(point_t a, point_t b, point_t c)
{
  double v1_x = a.x-b.x;
  double v1_y = a.y-b.y;
  double v2_x = c.x-b.x;
  double v2_y = c.y-b.y;
  double v1 = sqrt(v1_x*v1_x+v1_y*v1_y);
  double v2 = sqrt(v2_x*v2_x+v2_y*v2_y);
  if (v1>0 && v2>0)
    return (v1_x*v2_x+v1_y*v2_y)/(v1*v2);
  else
    return FLT_MAX;
}

// clockwise actual angle is considered positive
std::pair<double, double> find_rotation(point_t top, point_t left, point_t bottom, point_t right)
{
  double top_angle = fabs(cos_angle_between_points(left,top,right));
  double right_angle = fabs(cos_angle_between_points(top,right,bottom));
  double bottom_angle = fabs(cos_angle_between_points(right,bottom,left));
  double left_angle = fabs(cos_angle_between_points(bottom,left,top));
  double dx,dy;
  if (top_angle<right_angle && top_angle<left_angle && top_angle<bottom_angle)
    {
      // angle at the top point is the closest to 90 degrees
      if (top.x-left.x < right.x - top.x)
	{
	  // rotation is clockwise
	  dx = top.x-left.x;
	  dy = left.y-top.y;
	}
      else
	{
	  // rotation is counter-clockwise
	  dx = top.x-right.x;
	  dy = right.y-top.y;
	}
    }
  else if (left_angle<right_angle && left_angle<top_angle && left_angle<bottom_angle)
    {
      // angle at the left is the closest to 90 degrees
      if (bottom.y-left.y < left.y - top.y)
	{
	  // clockwise
	  dx = top.x-left.x;
	  dy = left.y-top.y;
	}
      else
	{
	  // counter-clockwise
	  dx = left.x-bottom.x;
	  dy = bottom.y-left.y;
	}
    }
  else if (right_angle<top_angle && right_angle<left_angle && right_angle<bottom_angle)
    {
      // angle at the right is the closest to 90 degrees
      if (right.y - top.y < bottom.y - right.y)
	{
	  // clockwise
	  dx = right.x-bottom.x;
	  dy = bottom.y-right.y;
	}
      else
	{
	  // counter-clockwise
	  dx = top.x-right.x;
	  dy = right.y-top.y;
	}
    }
  else
    {
      // assume that angle at the bottom is the closest to 90 degrees
      if (right.x-bottom.x < bottom.x-left.x)
	{
	  // clockwise
	  dx = right.x-bottom.x;
	  dy = bottom.y-right.y;
	}
      else 
	{
	  // counter-clockwise
	  dx = left.x - bottom.x;
	  dy = bottom.y - left.y;
	}
    }
  double s = dx / sqrt(dx*dx+dy*dy);
  double c = dy / sqrt(dx*dx+dy*dy);
  return(std::make_pair(s, c));
}

int border_count_in_rotated_frame(
    std::vector<std::vector<point_t> >::iterator m,
    point_t top_point, point_t left_point, point_t bottom_point,point_t right_point,
    std::pair <double, double> &sin_cos)
{
  double s = -sin_cos.first;
  double c = sin_cos.second;

  double x1 = top_point.x*c-top_point.y*s;
  double y1 = top_point.x*s+top_point.y*c;
  double x2 = left_point.x*c-left_point.y*s;
  double y2 = left_point.x*s+left_point.y*c;
  double x3 = bottom_point.x*c-bottom_point.y*s;
  double y3 = bottom_point.x*s+bottom_point.y*c;
  double x4 = right_point.x*c-right_point.y*s;
  double y4 = right_point.x*s+right_point.y*c;

  double left = std::min(std::min(x1, x2), std::min(x3, x4));
  double right = std::max(std::max(x1, x2), std::max(x3, x4));
  double top = std::min(std::min(y1, y2), std::min(y3, y4));
  double bottom = std::max(std::max(y1, y2), std::max(y3, y4));

  int border_count = 0;
  for (std::vector<point_t>::iterator p = m->begin(); p != m->end(); p++)
    {
      double x = c*(p->x) - s*(p->y);
      double y = s*(p->x) + c*(p->y);

      if ((x - left)<2 || (right - x) < 2 || (y - top) < 2 || (bottom - y) < 2)
	border_count++;
    }
  return border_count;
}


void remove_tables(std::vector<std::list<point_t> > &segments, std::vector<std::vector<point_t> > &margins,
                   unsigned int size)
{
  std::vector<std::list<point_t> >::iterator s;
  std::vector<std::vector<point_t> >::iterator m;
  s = segments.begin();
  m = margins.begin();

  while (s != segments.end() && m != margins.end())
    {
      if (m->size() <= size)
        {
          s++;
          m++;
          continue;
        }

      int top = INT_MAX, left = INT_MAX, bottom = 0, right = 0;
      point_t left_point,top_point,right_point,bottom_point;
      for (std::vector<point_t>::iterator p = m->begin(); p != m->end(); p++)
        {
          if (p->x < left)
	    {
	      left = p->x;
	      left_point = *p;
	    }
          if (p->x > right)
	    {
	      right = p->x;
	      right_point = *p;
	    }
          if (p->y < top)
	    {
	      top = p->y;
	      top_point = *p;
	    }
          if (p->y > bottom)
	    {
	      bottom = p->y;
	      bottom_point = *p;
	    }
        }

      double aspect = FLT_MAX;
      if (right != left)
        aspect = 1. * (bottom - top) / (right - left);
      if (aspect >= MAX_ASPECT || aspect <= 1./MAX_ASPECT)
        {
          s++;
          m++;
          continue;
        }

      double area = s->size();
      double square_area = (bottom - top+1) * (right - left+1);
      double ratio = 0;
      if (square_area != 0)
	ratio = area / square_area;
      if (ratio > MAX_RATIO || ratio == 0)
	{
          s++;
          m++;
          continue;
        }
      std::pair<double, double> zero_angle = std::make_pair(0, 1);
      int border_count = border_count_in_rotated_frame(m,top_point,left_point,bottom_point,right_point,zero_angle);

      if (PARTS_IN_MARGIN*border_count > BORDER_COUNT)
        {
          s = segments.erase(s);
          m = margins.erase(m);
        }
      else
        {
	  // perform rotation
          std::pair<double, double> sin_cos = find_rotation(top_point,left_point,bottom_point,right_point);
	  int rotated_border_count = border_count_in_rotated_frame(m,top_point,left_point,bottom_point,right_point,sin_cos);

	    if (PARTS_IN_MARGIN*rotated_border_count > BORDER_COUNT && fabs(sin_cos.first)<sin(10*PI/180) && sin_cos.second>cos(10*PI/180))
	    {
	      s = segments.erase(s);
	      m = margins.erase(m);
	    }
	  else
	    {
	      s++;
	      m++;
	    }
        }
    }
}

bool bulge(const point_t tail, const point_t head, const std::list<point_t> & seg)
{
  bool r = false;
  std::vector<int> y(std::max(abs(head.x - tail.x), abs(head.y - tail.y)) + 1, 0);
  int n=y.size();

  if (n<10) return false;

  bool horizontal = false;
  if (abs(head.x-tail.x)>abs(head.y-tail.y)) horizontal = true;

  for (std::list<point_t>::const_iterator p = seg.begin(); p != seg.end(); p++)
    {
      int d;
      if (horizontal)
	d=abs(p->x-tail.x);
      else
	d=abs(p->y-tail.y);
      if (d<n) y[d]++;
    }

  int top=0;
  int pos=0;
  for (int i=0; i<n; i++)
    if (y[i]>top)
      {
	top = y[i];
	pos = i;
      }

  /*for (int i=0; i<n; i++)
    std::cout<<y[i]<<" ";
    std::cout<<std::endl;*/
  
  if (pos<3) return false;
  int midpoint = std::min(int(0.75 * n), pos - 3);
  if (midpoint<0) return false;

  double avg=0;
  for (int i=0; i<midpoint; i++)
    avg +=y[i];
  avg /=int(midpoint);

  bool flat = true;
  for (int i=2; i<midpoint; i++)
    if (fabs(y[i]-avg)>2) flat = false;

  bool left = true;
  for (int i=pos-1; i>=midpoint; i--)
    if (y[i]>y[i+1]+2) 
      left=false;
  bool right = true;
  for (int i=pos+1; i<n; i++)
    if (y[i]>y[i-1]+2) right=false;
  bool peak = true;
  if (top<1.5*avg || top-avg<2 || n-pos<3 || top-y[n-1]<2 || pos<n/2 || pos<5) peak = false;

  //std::cout<<flat<<" "<<left<<" "<<right<<" "<<peak<<" "<<pos<<std::endl;

  return flat && left && right && peak;
}


void find_arrows_pluses(std::vector<std::vector<point_t> > &margins,
                        std::vector<std::list<point_t> > &segments,
                        std::vector<arrow_t> &arrows,
                        std::vector<plus_t> &pluses)
{
  const int len=50;
  for (int i=0; i<margins.size(); i++)
    {
      std::vector<int> hist(len,0);
      int top_pos=0;
      int top_value=0;
      point_t head, tail,center;
      int min_x, min_y, max_x, max_y;

      if (segments[i].size()>1000)
	build_hist(margins[i],hist,len,top_pos,top_value,head,tail,center,min_x, min_y, max_x, max_y);
      else
	build_hist(segments[i],hist,len,top_pos,top_value,head,tail,center,min_x, min_y, max_x, max_y);

      if (top_value>5)
	{
          std::vector<int> peaks(1,top_pos);
          std::vector<int> values(1,top_value);
	  for (int k=1; k<len;k++)
	    {
	      int pos=k+top_pos;
	      if (pos>=len) pos -= len;
	      int after=pos+1;
	      int before=pos-1;
	      int before2 = pos - 2;
	      if (after>=len) after -=len;
	      if (before<0) before +=len;
	      if (before2<0) before2 +=len;
	      if ( (hist[before]<hist[pos] || (hist[before] == hist[pos] && hist[before2] < hist[pos]))
		  && hist[after]<hist[pos] && hist[pos]>=top_value/2)  // find all peaks at least half as high as the top-most
		{
		  peaks.push_back(pos);
		  values.push_back(hist[pos]);
		}
	    }
	  if (peaks.size() == 2   && abs(len/2 - abs(peaks[1]-peaks[0]))<=1)  // only two peaks are present at 180 degrees
	    {
	      bool ba=bulge(tail,head,segments[i]);
	      bool bb=bulge(head,tail,segments[i]);
	      double l = distance(tail.x, tail.y, head.x, head.y);
	      //std::cout<<tail.x<<" "<<tail.y<<" "<<ba<<" " << bb << std::endl;
	      if ((ba || bb) && l > 2 * MAX_FONT_HEIGHT)
		{
		  // we found an arrow!
		  //std::cout<<tail.x<<" "<<tail.y<<std::endl;
		  arrow_t arrow(head,tail,min_x,min_y,max_x,max_y);
		  if (bb)
		    {
		      arrow.head = tail;
		      arrow.tail = head;
		    }
		  arrows.push_back(arrow);
		  margins[i].clear();
		  segments[i].clear();
		}
	    }
	  //std::cout << center.x << " " << center.y << std::endl;
	  //std::cout << "  " << values[0] << " " << values[1] << " " << values[2] <<" "<< values[3] << std::endl;
	  if (peaks.size() == 4  && (double(values[1])/values[0]>0.8 || values[0]-values[1]<=5)  && (double(values[2])/values[0]>0.8  || values[0]-values[2]<=5) && (double(values[3])/values[0]>0.8 || values[0]-values[3]<=5))
	    {
	      bool first=false, second=false, third=false, fourth=false;
	      for (int j=0; j<4; j++)
		{
		  if (peaks[j] <= 1 || peaks[j] >=len-2) first=true;
		  if (abs(len/4-peaks[j])<=1) second=true;
		  if (abs(len/2-peaks[j])<=1) third=true;
		  if (abs(3*len/4-peaks[j])<=1) fourth=true;
		}
	      /* for(int k=0; k<len; k++)
		{
		  std::cout << hist[k] << " ";
		}
		std::cout << std::endl;*/
	      for (int j=0; j<peaks.size(); j++)          // check outside of the peaks is essentially zero
		for (int k=peaks[j]-2; k<=peaks[j]+2; k++)
		  {
		    int kk=k;
		    if (kk<0) kk += len;
		    if (kk>=len) kk -=len;
		    hist[kk]=0;
		  }

	      bool low=true;
	      /* for(int k=0; k<len; k++)
		{
		  if (hist[k]>3) low=false;
		  }*/
	      //std::cout << "  " << first <<" "<<second << " " << third << " " << low << std::endl;
	      if (first && second && third && fourth && low)
		{
		  // we found a plus!
		  plus_t plus;
		  plus.center = center;
		  plus.min_x = min_x;
		  plus.min_y = min_y;
		  plus.max_x = max_x;
		  plus.max_y = max_y;
		  pluses.push_back(plus);
		}
	    }

	}
    }
  std::vector<std::vector<point_t> >::iterator k = margins.begin();
  while (k!=margins.end())
    {
      if (k->empty()) k = margins.erase(k);
      else k++;
    }
  std::vector<std::list<point_t> >::iterator l = segments.begin();
  while (l!=segments.end())
    {
      if (l->empty()) l = segments.erase(l);
      else l++;
    }

}

bool comp_labels(const label_t &left, const label_t &right)
{
  if (left.x2 < right.x1)
    return (true);
  if (std::max(left.y1, left.y2) < std::min(right.y1, right.y2))
    return (true);
  return (false);
}

int comp_labels_int(const void *l, const void *r)
{
  label_t *left = (label_t *) l;
  label_t *right = (label_t *) r;
  if (left->x2 < right->x1)
    return (-1);
  if (std::max(left->y1, left->y2) < std::min(right->y1, right->y2))
    return (-1);
  return (1);
}

std::string ocr_agent_strings(const std::vector<std::list<point_t> > &agents, const Image &image,
                              double threshold, const ColorGray &bgColor, bool verbose)
{
  std::string agent_string;
  std::vector<letters_t> letters;
  for (int i=0; i<agents.size(); i++)
    {
      int left=INT_MAX;
      int right=0;
      int top=INT_MAX;
      int bottom=0;
      for (std::list<point_t>::const_iterator a = agents[i].begin(); a != agents[i].end(); a++)
	{
	  if (a->x<left) left=a->x;
	  if (a->x>right) right=a->x;
	  if (a->y<top) top=a->y;
	  if (a->y>bottom) bottom=a->y;
	}
      if ((bottom - top) <= 2*MAX_FONT_HEIGHT && (right - left) <= 2*MAX_FONT_WIDTH && (bottom - top) > MIN_FONT_HEIGHT)
            {
              char label = 0;
              label = get_atom_label(image, bgColor, left, top, right, bottom, threshold, (right + left) / 2, top, true, verbose);
              if (label != 0)
                {
                  letters_t lt;
		  lt.a=label;
		  lt.x = (left + right) / 2;
		  lt.y  = (top + bottom) / 2;
		  lt.r  = distance(left, top, right, bottom) / 2;
		  lt.min_x = left;
		  lt.max_x = right;
		  lt.min_y = top;
		  lt.max_y = bottom;
		  lt.free = true;
                  letters.push_back(lt);
		}
	    }
    }
  std::vector<label_t> label;
  assemble_labels(letters, letters.size(), label);

  //sort(label.begin(),label.end(),comp_labels);

  if (!label.empty())
    qsort(&label[0],label.size(),sizeof(label_t),comp_labels_int);

  for (int i=0; i<label.size(); i++)
    if (!label[i].a.empty())
      agent_string += " "+label[i].a;
  return (agent_string);
}

void find_agent_strings(std::vector<std::vector<point_t> > &margins,
                        std::vector<std::list<point_t> > &segments, std::vector<arrow_t> &arrows,
			const Image &image, double threshold, const ColorGray &bgColor, bool verbose)
{
  for (int i=0; i<arrows.size(); i++)
    {
      std::vector<std::vector<point_t> > agent_margins;
      std::vector<std::list<point_t> > agents;

      double l=distance(arrows[i].tail.x,arrows[i].tail.y,arrows[i].head.x,arrows[i].head.y);
      bool found=false;
      for (int j=0; j<margins.size(); j++)
	{
	  bool close=false;
	  bool within=true;
	  for (int k=0; k<margins[j].size(); k++)
	    {
	      //     if (fabs(distance_from_bond_y(arrows[i].tail.x,arrows[i].tail.y,arrows[i].head.x,arrows[i].head.y,margins[j][k].x,margins[j][k].y))<MAX_FONT_HEIGHT) close=true;
	      //double d=distance_from_bond_x_a(arrows[i].tail.x,arrows[i].tail.y,arrows[i].head.x,arrows[i].head.y,margins[j][k].x,margins[j][k].y);
	      //if (d<0 || d>l) within=false;
	      if (fabs(distance((arrows[i].tail.x+arrows[i].head.x)/2,(arrows[i].tail.y+arrows[i].head.y)/2,margins[j][k].x,margins[j][k].y))<MAX_FONT_HEIGHT) close=true;
	    }
	  if (close && within)
	    {
	      agents.push_back(segments[j]);
	      agent_margins.push_back(margins[j]);
	      margins[j].clear();
	      segments[j].clear();
	      found=true;
	    }
	}

      while (found)
	{
	  found=false;

          std::vector<std::vector<point_t> >::iterator k = margins.begin();
	  while (k!=margins.end())
	    {
	      if (k->empty()) k = margins.erase(k);
	      else k++;
	    }
          std::vector<std::list<point_t> >::iterator l = segments.begin();
	  while (l!=segments.end())
	    {
	      if (l->empty()) l = segments.erase(l);
	      else l++;
	    }

	  for (int j=0; j<margins.size(); j++)
	    {
	      bool close=false;
	      for (int m=0; m<agent_margins.size(); m++)
		for (int k=0; k<margins[j].size(); k++)
		  for (int p=0; p<agent_margins[m].size(); p++)
		    if (distance(agent_margins[m][p].x,agent_margins[m][p].y,margins[j][k].x,margins[j][k].y)<MAX_FONT_HEIGHT/2) close=true;
	      if (close)
		{
		  agents.push_back(segments[j]);
		  agent_margins.push_back(margins[j]);
		  margins[j].clear();
		  segments[j].clear();
		  found=true;
		}
	    }
	}
      arrows[i].agent=ocr_agent_strings(agents,image,threshold,bgColor,verbose);
    }
}

std::vector<unsigned char> dilation(const std::vector<unsigned char> &image, int w, int h, unsigned int n)
{ 
  std::vector<unsigned char> out(w*h, 0);
  int before = (n - 1) / 2;
  int after = n - 1 - before;
  for(int i = 0; i < h; ++i)
    {
      int cache = 0;
      int j = 0;
      for (int ii = i - before; ii <= i + after; ++ii)
	for (int jj = j - before; jj <= j + after; ++jj)
	  if (ii >= 0 && ii < h && jj >= 0 && jj < w)
	    cache += image[ii*w + jj]; 
      if (cache > 0)
	out[i*w + j] = 1;
      ++j;
      for (; j < w; ++j)
	{
	  for (int ii = i - before; ii <= i + after; ++ii)
	    {
	      int jj = j - before; 
	      if (ii >= 0 && ii < h && jj >= 0 && jj < w)
		{
		  cache -= image[ii*w + jj]; 
		}	      
	      jj = j + after;
	      if (ii >= 0 && ii < h && jj >= 0 && jj < w)
		{
		  cache += image[ii*w + jj]; 
		}
	    }
           if (cache > 0)
	     out[i*w + j] = 1;
	}
    }
  return out;
}

std::vector<std::vector<unsigned char>> get_connected_components(const std::vector<unsigned char> &image, int w, int h)
{
  std::vector<bool> visited(w*h, false);
  std::vector<std::vector<unsigned char>> segments;
  for(int i = 0; i < h; ++i)
    for (int j = 0; j < w; ++j)
      {
	unsigned int pos = j + i * w;
	unsigned char c = image[pos];
	std::queue<unsigned int> bag;
	if (c == 1 && !visited[pos])
	  {
	    bag.push(pos);
	    visited[pos] = true;
	    std::vector<unsigned char> segment(w*h, 0);
	    while (!bag.empty())
	      {
		unsigned int current = bag.front();
		bag.pop();
		segment[current] = 1;
		int row = current / w;
		int col = current % w;
		for (int ii = row - 1; ii <= row + 1; ++ii)
		  for (int jj = col - 1; jj <= col + 1; ++jj)
		    if (ii >= 0 && ii < h && jj >= 0 && jj < w)
		      {
			pos = jj + ii * w;
			c = image[pos];
			if (c == 1 && !visited[pos])
			  {
			    bag.push(pos);
			    visited[pos] = true;
			  }
		      }
	      }
	    segments.push_back(segment);
	  }
      }
  return segments;
}

void save_connected_components(const std::vector<std::vector<unsigned char>> &images, const std::vector<std::list<point_t> > &segments,
			       int w, int h,
			       std::list<std::list<std::list<point_t> > > &explicit_clusters)
{
  explicit_clusters.clear();
  for (const auto &image : images)
    {
      std::list<std::list<point_t> > set_of_segments;
      int x1 = w;
      int y1 = h;
      int x2 = 0;
      int y2 = 0;
      int num_points = 0;
      for (const auto &segment : segments)
	{
	  point_t point = segment.front();
	  int pos = point.x + point.y * w;
	  if (image[pos] == 0)
	    continue;
	  for (auto pos : segment)
	    {
	      x1 = std::min(x1, pos.x);
	      y1 = std::min(y1, pos.y);
	      x2 = std::max(x2, pos.x);
	      y2 = std::max(y2, pos.y);
	    }
	  num_points += segment.size();
	  set_of_segments.push_back(segment);      
	}      
      if (set_of_segments.empty())
	continue;
      int new_w = x2 - x1 + 1;
      int new_h = y2 - y1 + 1;
      if (new_w < MAX_FONT_WIDTH || new_h < MAX_FONT_HEIGHT)
	continue;
      double aspect = static_cast<double>(std::max(new_w, new_h)) / std::min(new_w, new_h);
      double fill = static_cast<double>(num_points) / (new_w * new_h);
      if (aspect > MAX_ASPECT || fill > MAX_RATIO)
	continue;
      explicit_clusters.push_back(set_of_segments);
    }
}

void find_segments(
    const Image &image, double threshold, const ColorGray &bgColor, bool adaptive, bool is_reaction,
    std::vector<arrow_t> &arrows, std::vector<plus_t> &pluses, bool keep, bool verbose,
    std::list<std::list<std::list<point_t> > > &explicit_clusters, unsigned int segment_mask_size)
{
  explicit_clusters.clear();
  std::vector<std::list<point_t> > segments;
  std::vector<std::vector<point_t> > margins;
  

  // 1m34s

  find_connected_components(image, threshold, bgColor, segments, margins, adaptive);

  if (verbose)
    std::cout << "Number of segments: " << segments.size() << '.' << std::endl;

  if (segments.size() > MAX_SEGMENTS)
    {
      segments.clear();
      margins.clear();
    }
  if (is_reaction)
    {
      find_arrows_pluses(margins,segments,arrows, pluses);
      find_agent_strings(margins,segments,arrows,image,threshold,bgColor,verbose);
    }
  
  remove_separators(segments, margins, SEPARATOR_ASPECT, SEPARATOR_AREA);      
  remove_tables(segments, margins, SEPARATOR_AREA);
  // 2m22s
  
  if (keep)
    {
      std::list<int> new_cluster;
      for (unsigned int s = 0; s < margins.size(); s++)
	{
        new_cluster.push_back(s);       
      }
      std::list<std::list<int> > clusters;
      clusters.push_back(new_cluster);
      build_explicit_clusters(clusters, segments, explicit_clusters);
      int top = 0, left = INT_MAX, bottom = INT_MAX, right = 0;
      int num_points = 0;
      for (const auto &c : explicit_clusters)
	for (const auto &s : c)
	  for (const auto &p : s)
	    {
	      left = std::min(left, p.x);
	      right = std::max(right, p.x);
	      top = std::max(top, p.y);
	      bottom = std::min(bottom, p.y);
	      ++num_points;
	    }
      int new_w = right - left + 1;
      int new_h = top - bottom + 1;
      if (new_w < MAX_FONT_WIDTH || new_h < MAX_FONT_HEIGHT)
	{
	  explicit_clusters.clear();
	  return;
	}
      double aspect = static_cast<double>(std::max(new_w, new_h)) / std::min(new_w, new_h);
      double fill = static_cast<double>(num_points) / (new_w * new_h);
      if (aspect > MAX_ASPECT || fill > MAX_RATIO)
	{
	  explicit_clusters.clear();
	  return;
	} 
    }
  else
    {
      int w = image.columns();
      int h = image.rows();
      std::vector<unsigned char> in(h*w, 0);
      for (const auto &s : segments)
	for (const auto &p : s)
	  in[p.y * w + p.x] = 1;
      
      std::vector<unsigned char> out = dilation(in, w, h, segment_mask_size);

      /*std::vector<unsigned char> out1(out);
      for (auto &c : out1)
	c = 255 - 255*c;
      std::ofstream outfile("tmp.pgm", std::ios::binary );
      outfile << "P5" << std::endl;
      outfile << image.columns() << " " << image.rows() << std::endl;
      outfile << 255 << std::endl;
      std::copy(out1.begin(), out1.end(), std::ostreambuf_iterator<char>(outfile));
      */
      std::vector<std::vector<unsigned char>> split_images = get_connected_components(out, image.columns(),  image.rows());
      save_connected_components(split_images, segments, w, h, explicit_clusters);
    }
}

void find_box_size(const std::set<std::pair<int,int> > &set1, int &x1, int &x2, int &y1, int &y2)
{
  x1 = INT_MAX; y1 = INT_MAX; x2 = 0; y2 = 0;
  for (std::set<std::pair<int, int> >::const_iterator it = set1.begin(); it != set1.end(); ++it)
    {
      int px = it->first;
      int py = it->second;
      if (px < x1)
	x1 = px;
      if (px > x2)
	x2 = px;
      if (py < y1)
	y1 = py;
      if (py > y2)
	y2 = py;
    }
}

bool check_possible_bracket(std::set<std::pair<int, int> > &set1,
                            std::vector<std::vector<bool> > &global_pic,
                            int left, int right, int top, int bottom, int i)
{
  bool res = false;
  int x1,x2,y1,y2;
  find_box_size(set1,x1,x2,y1,y2);
  if (set1.size() <= 10 || (x2 - x1) <= 3 || (y2 - y1) <= 20 )
    return false;
  std::vector<std::vector<int> > tmp(x2 - x1 + 1, std::vector<int>(y2 - y1 + 1, 0));
  int startx = 0;
  int starty = INT_MAX;
  for (std::set<std::pair<int, int> >::const_iterator it = set1.begin(); it != set1.end(); ++it)
    {
      int px = it->first;
      int py = it->second;
      tmp[px-x1][py-y1] = 1;
      if (starty > py - y1)
	{
	  startx = px - x1;
	  starty = py - y1;
	}
    }
  int middle = (y2 - y1) / 2;
  std::vector<std::pair<int, int> > margin(1, std::make_pair(startx, starty));
  std::set<std::pair<int, int> > set2;
  while (!margin.empty())
    {
      startx = margin.back().first;
      starty = margin.back().second;
      margin.pop_back();
      int reflect = middle + (middle - starty);

      if (reflect >= 0 && reflect <= y2-y1 && tmp[startx][reflect] != 0)
	{
          set2.insert(std::make_pair(startx + x1, starty + y1));
	}
      tmp[startx][starty] = -1;
      for (int ii = startx-1; ii  <= startx+1; ii++)
	if (ii >= 0 && ii < tmp.size())
	  for (int j = starty-1; j <= starty+1; j++)
	    if (j >= 0 && j < tmp[ii].size() && tmp[ii][j] == 1)
	      {
		tmp[ii][j] = 2;
		margin.push_back(std::make_pair(ii, j));
	      }
    }
  swap(set1,set2);
  find_box_size(set1,x1,x2,y1,y2);
  if (set1.size() > 10 && (x2 - x1) > 3 && (y2 - y1) > 20 && (x2 - x1) < (y2 - y1))
    {
      int x = x2 - x1 + 1;
      int y = y2 - y1 + 1;
      int f = 1;
      // if (y > 40)
      //	f = y / 40;
      x /= f;
      y /= f;

      //      cout << x1 <<" " << y1 <<" "<<x2<<" "<<y2<<" " << i << endl;

      unsigned char *pic = (unsigned char *) malloc(x * y);
      for (int j = 0; j < x * y; j++)
	pic[j] = 255;
      for (std::set<std::pair<int, int> >::const_iterator it = set1.begin(); it != set1.end(); ++it)
	{
	  int px = it->first;
	  int py = it->second;
	  if ((py - y1) / f < y && (px - x1) / f < x && (py - y1) % f == 0 && (px - x1) % f == 0)
	    pic[((py - y1) / f) * x + (px - x1) / f] = 0;
	}

      res = detect_bracket(x, y, pic);
    }
  return res;
}

void remove_brackets(int left, int right, int top, int bottom,
                     std::list<std::list<std::list<point_t> > >::iterator c,
                     std::set<std::pair<int, int> > &brackets)
{
  std::vector<std::vector<bool> > tmp(right - left + 1, std::vector<bool> (bottom - top + 1, false));
  std::vector<std::vector<bool> > global_pic(right - left + 1, std::vector<bool> (bottom - top + 1, false));

  for (std::list<std::list<point_t> >::const_iterator s = c->begin(); s != c->end(); s++)
    for (std::list<point_t>::const_iterator p = s->begin(); p != s->end(); p++)
      global_pic[p->x - left][p->y - top] = true;

  //Image t(Geometry(right + 1, bottom + 1), "white");

 for (int i = left + FRAME; i < right - FRAME; i++)
    {
      for (std::list<std::list<point_t> >::const_iterator s = c->begin(); s != c->end(); s++)
	{
          std::set<std::pair<int, int> > set1;
	  for (std::list<point_t>::const_iterator p = s->begin(); p != s->end(); p++)
	    {
	      if ( i + (i - p->x) < right &&
		   i + (i - p->x) - left - 1 >= 0 &&  p->y - top - 1 >= 0 &&
		   i + (i - p->x) - left + 1 < right - left + 1 && p->y - top + 1 < bottom - top + 1 &&
		   global_pic[p->x - left][p->y - top] &&
		   (global_pic[i + (i - p->x) - left][p->y - top] ||
		    global_pic[i + (i - p->x) - left + 1][p->y - top] ||
		    global_pic[i + (i - p->x) - left - 1][p->y - top] ||
		    global_pic[i + (i - p->x) - left][p->y - top + 1] ||
		    global_pic[i + (i - p->x) - left][p->y - top - 1] ||
		    global_pic[i + (i - p->x) - left - 1][p->y - top - 1] ||
		    global_pic[i + (i - p->x) - left - 1][p->y - top + 1] ||
		    global_pic[i + (i - p->x) - left + 1][p->y - top - 1] ||
		    global_pic[i + (i - p->x) - left + 1][p->y - top + 1]  )
		   && (p->x < i - 40) )
		{
                  set1.insert(std::make_pair(p->x, p->y));
		}
	    }

	  if (check_possible_bracket(set1, global_pic, left, right, top, bottom, i))
	    {
              for (std::set<std::pair<int, int> >::const_iterator it = set1.begin(); it != set1.end(); ++it)
		{
		  int ox = it->first;
		  int oy = it->second;
		  int px = i + (i - ox);
		  brackets.insert(*it);
		  brackets.insert(std::make_pair(px, oy));
		  //t.pixelColor(ox, oy, "black");
		  //t.pixelColor(px, oy, "black");
		}
	    }
	}
    }
 //t.write("t.png");
}


int prune_clusters(std::list<std::list<std::list<point_t> > > &clusters, std::vector<box_t> &boxes,
                   std::set<std::pair<int, int> > &brackets)
{
  int n_boxes = 0;
  std::list<std::list<std::list<point_t> > >::iterator c = clusters.begin();

  while (c != clusters.end())
    {
      int top = INT_MAX, left = INT_MAX, bottom = 0, right = 0;
      for (std::list<std::list<point_t> >::const_iterator s = c->begin(); s != c->end(); s++)
        {
          for (std::list<point_t>::const_iterator p = s->begin(); p != s->end(); p++)
            {
	      left = std::min(left, p->x);
	      right = std::max(right, p->x);
	      top = std::min(top, p->y);
	      bottom = std::max(bottom, p->y);
            }
        }

      box_t b1;
      boxes.push_back(b1);
      boxes[n_boxes].x1 = left;
      boxes[n_boxes].y1 = top;
      boxes[n_boxes].x2 = right;
      boxes[n_boxes].y2 = bottom;
      
      remove_brackets(left, right, top, bottom, c, brackets);
      
      for (std::list<std::list<point_t> >::const_iterator s = c->begin(); s != c->end(); s++)
	for (std::list<point_t>::const_iterator p = s->begin(); p != s->end(); p++)
	  boxes[n_boxes].c.push_back(*p);
      c++;
      n_boxes++;
    }
  return (n_boxes);
}

