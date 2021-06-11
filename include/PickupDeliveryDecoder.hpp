#ifndef LAB_MC658_PICKUPDELIVERYDECODER_HPP
#define LAB_MC658_PICKUPDELIVERYDECODER_HPP

#include "pickup_delivery_utils.hpp"
#include <algorithm>
#include <list>
#include <vector>

class PickupDeliveryDecoder {
public:
  PickupDeliveryDecoder();
  ~PickupDeliveryDecoder();

  double decode(const std::vector<double> &chromosome) const;
  void decode(const std::vector<double> &chromosome, DNodeVector &Sol) const;
};

#endif // LAB_MC658_PICKUPDELIVERYDECODER_HPP
