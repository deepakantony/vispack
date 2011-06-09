// ******************************************************************
// **********************  SystemFitModel *******************
// ******************************************************************

template <class TProjection, class TModel>
float SystemFitModel<TProjection,TModel>::errorSq() const
{
  int n = numPoints();
  int i;
  float ret = 0.0f;
  for (i = 0; i < n; i++)
    ret += power(errorX(i), 2);
  return(ret);
};

template <class TProjection, class TModel>
VISVector SystemFitModel<TProjection,TModel>::dEdParams(int i) const
{
  cout << "SystemProjectionSquare::dEdParams(int i)" << endl;
  VISMatrix los = _projection.lineOfSight(i);
  cout << "los" << los << endl;
  VISArray< VISVector > intersects = _model.intersections(los);
  cout << "intersects.n " << intersects.n() << endl;
  int n = intersects.n();
  VISVector ret(_model.numParams()), point, normal;
  ret = 0.0f;
  for (int i = 0; i < n; i++)
    {
      point = intersects.peek(i);
      normal = _model.normal(point);
      cout << "point" << point << endl;
      cout << "normal" << normal << endl;
      cout << "dxdParams" << _model.dXdParams(point) << endl;
      ret += ((_model.dXdParams(point)).t())*normal;
    }
  return(ret);
};



//  // ******************************************************************
//  // **********************  SystemModel <templated>*******************
//  // ******************************************************************

//  template <class TProjection, class TModel>
//  float SystemModel::errorSq() const
//  {
//    int n = numPoints();
//    int i;
//    float ret = 0.0f;
//    for (i = 0; i < n; i++)
//      ret += power(errorX(i), 2);
//    return(ret);
//  }

//  VISVector SystemProjectionSquare::dEdParams(int i) const
//  {
//    cout << "SystemProjectionSquare::dEdParams(int i)" << endl;
//    VISMatrix los = _projection.lineOfSight(i);
//    cout << "los" << los << endl;
//    VISArray< VISVector > intersects = _square.intersections(los);
//    cout << "intersects.n " << intersects.n() << endl;
//    int n = intersects.n();
//    VISVector ret(_square.numParams()), point, normal;
//    ret = 0.0f;
//    for (int i = 0; i < n; i++)
//      {
//        point = intersects.peek(i);
//        normal = _square.normal(point);
//        cout << "point" << point << endl;
//        cout << "normal" << normal << endl;
//        cout << "dxdParams" << _square.dXdParams(point) << endl;
//        ret += ((_square.dXdParams(point)).t())*normal;
//      }
//    return(ret);
//  }
